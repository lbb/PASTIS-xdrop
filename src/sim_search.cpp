/**
 * @file
 *	sim_search.cpp
 *
 * @author
 *	Oguz Selvitopi
 *
 * @date
 *
 * @brief
 *	PASTIS pipeline for computing seq similarity matrix
 *
 * @todo
 *
 * @note
 *	
 */

#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../inc/align/IPumaAligner.hpp"
#include "../inc/align/AdeptBSWAligner.hpp"
#include "../inc/align/PWAlign.hpp"
#include "../inc/align/SeqanFullAligner.hpp"
#include "../inc/align/SeqanXdropAligner.hpp"
#include "../inc/Alphabet.hpp"
#include "../inc/DistFastaData.hpp"
#include "../inc/DistPWRunner.hpp"
#include "../inc/kmer/Kmer.hpp"
#include "../inc/kmer/KmerOps.hpp"
#include "../inc/ParallelOps.hpp"
#include "../inc/sim_search.hpp"
#include "../inc/sr.hpp"
#include "../inc/Types.hpp"
#include "../inc/util.hpp"



using std::string;	using std::shared_ptr;	using std::vector;	using std::move;
using std::pair;	using std::tuple;	using std::to_string;
using std::unordered_set;

extern shared_ptr<pastis::ParallelOps> parops; // global var




namespace
pastis
{




// computes all blocks and prunes based on even or odd index ids
template <typename T_OUT>
void
mult_aln_bl_idx
(
 	params_t						&params,
   	shared_ptr<DistFastaData>		 dfd,
 	PSpMat<MatrixEntry>::MPI_DCCols &A,
	PSpMat<MatrixEntry>::MPI_DCCols &AT
)
{
	parops->logger->log("Computing SpGEMM and alignment (blocked version).");
	
	// types
	using SR = typename std::conditional<
		std::is_same<T_OUT, CommonKmerLight>::value,
		KmerIntersectLight<MatrixEntry, T_OUT>,
		KmerIntersectLoc<MatrixEntry, T_OUT>
		>::type;
	using T_SpMat = typename PSpMat<T_OUT>::MPI_DCCols;
	using DER_IN = typename combblas::SpDCCols<uint64_t, MatrixEntry>;


	parops->tp->start_timer("seq-comm-wait");
	
	// wait till receive block seqs
	if (!dfd->is_ready())
	  	dfd->bl_wait();

	parops->tp->stop_timer("seq-comm-wait");



	parops->tp->start_timer("sparse");
	parops->tp->start_timer("block-split");
	
	// block spgemm instance (splits matrices)
	combblas::BlockSpGEMM<uint64_t, MatrixEntry, DER_IN,
						  MatrixEntry, DER_IN>
		bspgemm(A, AT, params.br, params.bc);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("block-split");


	// block matrix and combblas matrix offsets
	vector<pair<uint64_t, uint64_t>> cb_roffsets =
		cbmat_offsets(A.getnrow(), params.br, false);
	vector<pair<uint64_t, uint64_t>> cb_coffsets =
		cbmat_offsets(AT.getncol(), params.bc, true);


	// final similarity matrix els
	vector<uint64_t> rids;
	vector<uint64_t> cids;
	vector<float> vals;

	
	// aligner
	PWAlign *pwa = nullptr;
	if (params.pw_aln != params_t::PwAln::ALN_NONE)
	{
		if (params.pw_aln == params_t::PwAln::ALN_SEQAN_FULL)
			pwa = new SeqanFullAligner(params.gap_open, params.gap_ext);
		else if (params.pw_aln == params_t::PwAln::ALN_SEQAN_XDROP)
			pwa = new SeqanXdropAligner
				(params.gap_open, params.gap_ext, params.klength,
				 params.aln_seqan_xdrop, params.seed_count);
#ifdef ADEPT_ENABLED
		else if (params.pw_aln == params_t::PwAln::ALN_ADEPT_GPUBSW)
			pwa = new AdeptBSWAligner(params.gap_open, params.gap_ext);
#endif
		else if (params.pw_aln == params_t::PwAln::ALN_IPUMA)
			pwa = new IPumaAligner(params.gap_open, params.gap_ext, params.klength, params.seed_count);

		parops->tp->start_timer("construct-seqs");
		
		// raw fasta data to aligner-specific representation
		pwa->construct_seqs_bl(dfd);

		parops->tp->stop_timer("construct-seqs");
	}


	parops->tp->start_timer("mult-align");


	// blocked multiplication and alignment
	parops->tp->start_timer("sim:block_all");
	uint64_t tot_cnnz = 0;
	for (int bri = 0; bri<params.br; ++bri)
	{
		for (int bci = 0; bci<params.bc; ++bci)
		{
			parops->info("Block " + to_string(bri) + " " +
						 to_string(bci) + " begin");
			auto t_block = std::chrono::system_clock::now();
			
			parops->tp->start_timer("sparse");
			parops->tp->start_timer("mat-mult");

			auto t_begin = std::chrono::system_clock::now();
			
			uint64_t roffset, coffset;
			auto C =
				bspgemm.getBlockId<SR, T_OUT,
								   combblas::SpDCCols<uint64_t, T_OUT>>
								   (bri, bci, roffset, coffset);

			double t_block_pure =
				ms_t(std::chrono::system_clock::now()-t_begin).count();
			parops->info("Pure block took " +
						 to_string(t_block_pure/1e3) + " seconds");

			parops->info("Matrix mult complete");

			parops->tp->stop_timer("mat-mult");
			parops->tp->start_timer("prune");
			
			uint64_t cnnz = C.getnnz();
			parops->tp->mat_stats["c_nnz"] += cnnz;

			auto f_eo_idx =
				[&roffset, &coffset]
				(const tuple<uint64_t, uint64_t, T_OUT> &t)
				{
					uint64_t r = std::get<0>(t)+roffset;
					uint64_t c = std::get<1>(t)+coffset;
					bool is_r_odd = r & 0x1;
					bool is_c_odd = c & 0x1;
					return (r == c) ||
						(r > c ?
						 (is_r_odd ^ is_c_odd) :
						 (!(is_r_odd ^ is_c_odd)));
						 
				};
			C.PruneI(f_eo_idx);

			parops->tp->stop_timer("prune");
			parops->tp->stop_timer("sparse");

			parops->info("Symmetricity pruning complete");
			
			uint64_t cnnz_partial = C.getnnz();
			parops->tp->mat_stats["prune_sym"] += cnnz_partial;

			if (pwa)
				pw_aln_batch(dfd, C, *pwa, params, bri, bci);

			uint64_t cnnz_aln = C.getnnz();
			parops->info
				("Block " + to_string(bri) + " " + to_string(bci) + 
				 " computed C nnz " + to_string(cnnz) +
				 " before aln C (partial block) " + to_string(cnnz_partial) +
				 " after aln " + to_string(C.getnnz()) +
				 " elim " + 
				 to_string((static_cast<double>(cnnz-cnnz_aln)/
							cnnz)*100.0) + "%");
			tot_cnnz += cnnz_aln;
			parops->tp->mat_stats["sim_thr"] += cnnz_aln;

			if (C.seqptr()->getnnz() == 0)
				continue;

			parops->tp->start_timer("mat-formation");
			
			auto dcsc = C.seqptr()->GetDCSC();
			for (uint64_t i = 0; i < dcsc->nzc; ++i)
			{
				for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
				{
					rids.push_back(dcsc->ir[j]+cb_roffsets[bri].first);
					cids.push_back(dcsc->jc[i]+cb_coffsets[bci].first);
					vals.push_back(dcsc->numx[j].score_aln);
				}
			}

			parops->tp->stop_timer("mat-formation");

			#if PASTIS_DBG_LVL > 0
			parops->bytes_alloc += rids.capacity() * sizeof(uint64_t);
			parops->bytes_alloc += cids.capacity() * sizeof(uint64_t);
			parops->bytes_alloc += vals.capacity() * sizeof(float);
			parops->logger->log("approximate memory in usage " +
								gb_str(parops->bytes_alloc));
			#endif


			double t_block_cur =
				ms_t(std::chrono::system_clock::now()-t_block).count();
			parops->info
				("Block " + to_string(bri) + " " + to_string(bci) +
				 " took " + to_string(t_block_cur/1e3) + " seconds");
		}
	}

	parops->tp->stop_timer("mult-align");
	parops->tp->start_timer("sparse");
	parops->tp->start_timer("sim-mat-assemble");

	combblas::FullyDistVec<uint64_t, uint64_t> drows(rids, parops->grid);
	combblas::FullyDistVec<uint64_t, uint64_t> dcols(cids, parops->grid);
	combblas::FullyDistVec<uint64_t, float> dvals(vals, parops->grid);
	combblas::SpParMat<uint64_t, float,
					   combblas::SpDCCols<uint64_t, float>>
		S(A.getnrow(), AT.getncol(), drows, dcols, dvals);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("sim-mat-assemble");
	parops->tp->start_timer("sim-mat-io");

	if (!params.align_file.empty())
		S.ParallelWriteMM(params.align_file, true);

	parops->tp->stop_timer("sim-mat-io");
	
	parops->info("tot nnz in similarity mat " + to_string(S.getnnz()));


	parops->tp->stop_timer("main:sim->mult_align");


	if (pwa)
		delete pwa;
}



// computes all blocks and prunes based on even or odd index ids
// overlap version
template <typename T_OUT>
void
mult_aln_bl_idx_ovlp
(
 	params_t						&params,
   	shared_ptr<DistFastaData>		 dfd,
 	PSpMat<MatrixEntry>::MPI_DCCols &A,
	PSpMat<MatrixEntry>::MPI_DCCols &AT
)
{
	parops->logger->log("Computing SpGEMM and alignment (blocked version).");
	
	// types
	using SR = typename std::conditional<
		std::is_same<T_OUT, CommonKmerLight>::value,
		KmerIntersectLight<MatrixEntry, T_OUT>,
		KmerIntersectLoc<MatrixEntry, T_OUT>
		>::type;
	using T_SpMat = typename PSpMat<T_OUT>::MPI_DCCols;
	using DER_IN = typename combblas::SpDCCols<uint64_t, MatrixEntry>;

	parops->tp->start_timer("seq-comm-wait");

	// wait till receive block seqs
	if (!dfd->is_ready())
	  	dfd->bl_wait();

	parops->tp->stop_timer("seq-comm-wait");

	parops->tp->start_timer("sparse");
	parops->tp->start_timer("block-split");

	// block spgemm instance (splits matrices)
	combblas::BlockSpGEMM<uint64_t, MatrixEntry, DER_IN,
						  MatrixEntry, DER_IN>
		bspgemm(A, AT, params.br, params.bc);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("block-split");

	// block matrix and combblas matrix offsets
	vector<pair<uint64_t, uint64_t>> cb_roffsets =
		cbmat_offsets(A.getnrow(), params.br, false);
	vector<pair<uint64_t, uint64_t>> cb_coffsets =
		cbmat_offsets(AT.getncol(), params.bc, true);

	// final similarity matrix els
	vector<uint64_t> rids;
	vector<uint64_t> cids;
	vector<float> vals;
	
	// aligner
	PWAlign *pwa = nullptr;
	if (params.pw_aln != params_t::PwAln::ALN_NONE)
	{
		if (params.pw_aln == params_t::PwAln::ALN_SEQAN_FULL)
			pwa = new SeqanFullAligner(params.gap_open, params.gap_ext);
		else if (params.pw_aln == params_t::PwAln::ALN_SEQAN_XDROP)
			pwa = new SeqanXdropAligner
				(params.gap_open, params.gap_ext, params.klength,
				 params.aln_seqan_xdrop, params.seed_count);
#ifdef ADEPT_ENABLED
		else if (params.pw_aln == params_t::PwAln::ALN_ADEPT_GPUBSW)
			pwa = new AdeptBSWAligner(params.gap_open, params.gap_ext);
#endif
		else if (params.pw_aln == params_t::PwAln::ALN_IPUMA)
			pwa = new IPumaAligner(params.gap_open, params.gap_ext, params.klength, params.seed_count);

		parops->tp->start_timer("construct-seqs");

		// raw fasta data to aligner-specific representation
		pwa->construct_seqs_bl(dfd);

		parops->tp->stop_timer("construct-seqs");
	}


	// setup thread data for overlapping
	MultData<MatrixEntry, T_OUT> md(&bspgemm, &params);

	parops->tp->start_timer("mult-align");

	// blocked multiplication and alignment
	uint64_t tot_cnnz = 0, cur_block = 0;
	while (cur_block++ < params.br * params.bc)
	{
		int bri = md.bri_, bci = md.bci_;
		if (bci == -1)
			bci = 0;		
		parops->info("Block " + to_string(bri) + " " +
						 to_string(bci) + " begin");

		auto t_block = std::chrono::system_clock::now();

		// first but the other blocks will be overlapped with alignments
		if (md.bri_ == 0 && md.bci_ == -1)
			md.computeNext();

		parops->info("Matrix mult complete");

		parops->tp->start_timer("sparse");
		parops->tp->start_timer("prune");				

		auto &C = md.C_[md.curC_];

		uint64_t cnnz = C.getnnz();
		parops->tp->mat_stats["c_nnz"] += cnnz;

		uint64_t roffset = md.roffset_, coffset = md.coffset_;
		auto f_eo_idx =
			[&roffset, &coffset]
			(const tuple<uint64_t, uint64_t, T_OUT> &t)
			{
				uint64_t r = std::get<0>(t)+roffset;
				uint64_t c = std::get<1>(t)+coffset;
				bool is_r_odd = r & 0x1;
				bool is_c_odd = c & 0x1;
				return (r == c) ||
					(r > c ?
					 (is_r_odd ^ is_c_odd) :
					 (!(is_r_odd ^ is_c_odd)));

			};
		C.PruneI(f_eo_idx);
		
		parops->info("Symmetricity pruning complete");
		parops->tp->stop_timer("prune");
		parops->tp->stop_timer("sparse");		

		uint64_t cnnz_partial = C.getnnz();
		parops->tp->mat_stats["prune_sym"] += cnnz_partial;

		if (pwa)
			pw_aln_batch_ovlp(dfd, C, *pwa, params, md, bri, bci);

		uint64_t cnnz_aln = C.getnnz();
		parops->info
			("Block " + to_string(bri) + " " + to_string(bci) + 
			 " computed C nnz " + to_string(cnnz) +
			 " before aln C (partial block) " + to_string(cnnz_partial) +
			 " after aln " + to_string(C.getnnz()) +
			 " elim " + 
			 to_string((static_cast<double>(cnnz-cnnz_aln)/
						cnnz)*100.0) + "%");
		tot_cnnz += cnnz_aln;
		parops->tp->mat_stats["sim_thr"] += cnnz_aln;

		if (C.seqptr()->getnnz() == 0)
			continue;

		parops->tp->start_timer("mat-formation");

		auto dcsc = C.seqptr()->GetDCSC();
		for (uint64_t i = 0; i < dcsc->nzc; ++i)
		{
			for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
			{
				rids.push_back(dcsc->ir[j]+cb_roffsets[bri].first);
				cids.push_back(dcsc->jc[i]+cb_coffsets[bci].first);
				vals.push_back(dcsc->numx[j].score_aln);
			}
		}

		parops->tp->stop_timer("mat-formation");

		#if PASTIS_DBG_LVL > 0
		parops->bytes_alloc += rids.capacity() * sizeof(uint64_t);
		parops->bytes_alloc += cids.capacity() * sizeof(uint64_t);
		parops->bytes_alloc += vals.capacity() * sizeof(float);
		parops->logger->log("approximate memory in usage " +
							gb_str(parops->bytes_alloc));
		#endif


		double t_block_cur =
			ms_t(std::chrono::system_clock::now()-t_block).count();
		parops->info
			("Block " + to_string(bri) + " " + to_string(bci) +
			 " took " + to_string(t_block_cur/1e3) + " seconds");
	}

	parops->tp->stop_timer("mult-align");

	// parops->tp->stop_timer("total");
	// string stats = parops->tp->to_string();
	// if (parops->g_rank == 0)
	// 	std::cout << stats << endl;
	// parops->tp->start_timer("total");
	
	parops->tp->start_timer("sparse");
	parops->tp->start_timer("sim-mat-assemble");

	combblas::FullyDistVec<uint64_t, uint64_t> drows(rids, parops->grid);
	combblas::FullyDistVec<uint64_t, uint64_t> dcols(cids, parops->grid);
	combblas::FullyDistVec<uint64_t, float> dvals(vals, parops->grid);
	combblas::SpParMat<uint64_t, float,
					   combblas::SpDCCols<uint64_t, float>>
		S(A.getnrow(), AT.getncol(), drows, dcols, dvals);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("sim-mat-assemble");
	parops->tp->start_timer("sim-mat-io");

	if (!params.align_file.empty())
		S.ParallelWriteMM(params.align_file, true);

	parops->tp->stop_timer("sim-mat-io");
	
	parops->info("tot nnz in similarity mat " + to_string(S.getnnz()));


	if (pwa)
		delete pwa;
}



template <typename T_OUT>
void
mult_aln_bl_trg
(
 	params_t						&params,
   	shared_ptr<DistFastaData>		 dfd,
 	PSpMat<MatrixEntry>::MPI_DCCols &A,
	PSpMat<MatrixEntry>::MPI_DCCols &AT
)
{
	parops->logger->log("Computing SpGEMM and alignment (blocked version).");
	
	// types
	using SR = typename std::conditional<
		std::is_same<T_OUT, CommonKmerLight>::value,
		KmerIntersectLight<MatrixEntry, T_OUT>,
		KmerIntersectLoc<MatrixEntry, T_OUT>
		>::type;
	using T_SpMat = typename PSpMat<T_OUT>::MPI_DCCols;
	using DER_IN = typename combblas::SpDCCols<uint64_t, MatrixEntry>;


	parops->tp->start_timer("seq-comm-wait");
	
	// wait till receive block seqs
	if (!dfd->is_ready())
	  	dfd->bl_wait();
	
	parops->tp->stop_timer("seq-comm-wait");

	
	parops->tp->start_timer("sparse");
	parops->tp->start_timer("block-split");
	
	// block spgemm instance (splits matrices)
	combblas::BlockSpGEMM<uint64_t, MatrixEntry, DER_IN,
						  MatrixEntry, DER_IN>
		bspgemm(A, AT, params.br, params.bc);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("block-split");


	// block matrix and combblas matrix offsets
	vector<uint64_t> roffsets = move(bspgemm.getBlockOffsets(true));
	vector<uint64_t> coffsets = move(bspgemm.getBlockOffsets(false));
	vector<pair<uint64_t, uint64_t>> cb_roffsets =
		cbmat_offsets(A.getnrow(), params.br, false);
	vector<pair<uint64_t, uint64_t>> cb_coffsets =
		cbmat_offsets(AT.getncol(), params.bc, true);


	// final similarity matrix els
	vector<uint64_t> rids;
	vector<uint64_t> cids;
	vector<float> vals;

	
	// aligner
	PWAlign *pwa = nullptr;
	if (params.pw_aln != params_t::PwAln::ALN_NONE)
	{
		if (params.pw_aln == params_t::PwAln::ALN_SEQAN_FULL)
			pwa = new SeqanFullAligner(params.gap_open, params.gap_ext);
		else if (params.pw_aln == params_t::PwAln::ALN_SEQAN_XDROP)
			pwa = new SeqanXdropAligner
				(params.gap_open, params.gap_ext, params.klength,
				 params.aln_seqan_xdrop, params.seed_count);			
#ifdef ADEPT_ENABLED
	    	else if (params.pw_aln == params_t::PwAln::ALN_ADEPT_GPUBSW)
			pwa = new AdeptBSWAligner(params.gap_open, params.gap_ext);
#endif
		else if (params.pw_aln == params_t::PwAln::ALN_IPUMA)
			pwa = new IPumaAligner(params.gap_open, params.gap_ext, params.klength, params.seed_count);

				
		parops->tp->start_timer("construct-seqs");
		
		// raw fasta data to aligner-specific representation
		pwa->construct_seqs_bl(dfd);
		
		parops->tp->stop_timer("construct-seqs");		
	}


	parops->tp->start_timer("mult-align");

	// blocked multiplication and alignment
	uint64_t tot_cnnz = 0;
	for (int bri = 0; bri<params.br; ++bri)
	{
		for (int bci = 0; bci<params.bc; ++bci)
		{
			parops->info("Block " + to_string(bri) + " " +
						 to_string(bci) + " begin");
			auto t_block = std::chrono::system_clock::now();
			
			bool full_block = coffsets[bci] > (roffsets[bri+1]-1);
			// diag els skipped
			bool skip_block = roffsets[bri] >= (coffsets[bci+1]-1);
			bool partial_block = !full_block && !skip_block;

			if (skip_block)
				continue;

			parops->tp->start_timer("sparse");
			parops->tp->start_timer("mat-mult");

			auto t_begin = std::chrono::system_clock::now();

			uint64_t roffset, coffset;
			auto C =
				bspgemm.getBlockId<SR, T_OUT,
								   combblas::SpDCCols<uint64_t, T_OUT>>
								   (bri, bci, roffset, coffset);

			double t_block_pure =
				ms_t(std::chrono::system_clock::now()-t_begin).count();
			parops->info("Pure block took " +
						 to_string(t_block_pure/1e3) + " seconds");

			parops->tp->stop_timer("mat-mult");
			parops->tp->start_timer("prune");
			
			uint64_t cnnz = C.getnnz();
			parops->tp->mat_stats["c_nnz"] += cnnz;
			if (partial_block)	// get rid of lower triangular matrix elems
			{
				auto f_strict_upper =
					[&roffset, &coffset]
					(const tuple<uint64_t, uint64_t, T_OUT> &t)
					{
						return std::get<0>(t)+roffset >= std::get<1>(t)+coffset;
					};
				C.PruneI(f_strict_upper);
			}

			parops->tp->stop_timer("prune");
			parops->tp->stop_timer("sparse");

			uint64_t cnnz_partial = C.getnnz();
			parops->tp->mat_stats["prune_sym"] += cnnz_partial;

			if (pwa)
				pw_aln_batch(dfd, C, *pwa, params, bri, bci);

			uint64_t cnnz_aln = C.getnnz();
			parops->info
				("Block " + to_string(bri) + " " + to_string(bci) + 
				 " computed C nnz " + to_string(cnnz) +
				 " before aln C (partial block) " + to_string(cnnz_partial) +
				 " after aln " + to_string(C.getnnz()) +
				 " elim " + 
				 to_string((static_cast<double>(cnnz-cnnz_aln)/
							cnnz)*100.0) + "%");
			tot_cnnz += cnnz_aln;
			parops->tp->mat_stats["sim_thr"] += cnnz_aln;

			if (C.seqptr()->getnnz() == 0)
				continue;

			parops->tp->start_timer("mat-formation");
			
			auto dcsc = C.seqptr()->GetDCSC();
			for (uint64_t i = 0; i < dcsc->nzc; ++i)
			{
				for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
				{
					rids.push_back(dcsc->ir[j]+cb_roffsets[bri].first);
					cids.push_back(dcsc->jc[i]+cb_coffsets[bci].first);
					vals.push_back(dcsc->numx[j].score_aln);
				}
			}

			parops->tp->stop_timer("mat-formation");

			#if PASTIS_DBG_LVL > 0
			parops->bytes_alloc += rids.capacity() * sizeof(uint64_t);
			parops->bytes_alloc += cids.capacity() * sizeof(uint64_t);
			parops->bytes_alloc += vals.capacity() * sizeof(float);
			parops->logger->log("approximate memory in usage " +
								gb_str(parops->bytes_alloc));
			#endif

			double t_block_cur =
				ms_t(std::chrono::system_clock::now()-t_block).count();
			parops->info
				("Block " + to_string(bri) + " " + to_string(bci) +
				 " took " + to_string(t_block_cur/1e3) + " seconds");
		}
	}

	parops->tp->stop_timer("mult-align");
	parops->tp->start_timer("sparse");
	parops->tp->start_timer("sim-mat-assemble");

	combblas::FullyDistVec<uint64_t, uint64_t> drows(rids, parops->grid);
	combblas::FullyDistVec<uint64_t, uint64_t> dcols(cids, parops->grid);
	combblas::FullyDistVec<uint64_t, float> dvals(vals, parops->grid);
	combblas::SpParMat<uint64_t, float,
					   combblas::SpDCCols<uint64_t, float>>
		S(A.getnrow(), AT.getncol(), drows, dcols, dvals);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("sim-mat-assemble");
	parops->tp->start_timer("sim-mat-io");
	
	if (!params.align_file.empty())
		S.ParallelWriteMM(params.align_file, true);

	parops->tp->stop_timer("sim-mat-io");
	
	parops->info("tot nnz in similarity mat " + to_string(S.getnnz()));


	if (pwa)
		delete pwa;
}




template <typename T_OUT>
void
mult_aln_bl_trg_ovlp
(
 	params_t						&params,
   	shared_ptr<DistFastaData>		 dfd,
 	PSpMat<MatrixEntry>::MPI_DCCols &A,
	PSpMat<MatrixEntry>::MPI_DCCols &AT
)
{
	parops->logger->log("Computing SpGEMM and alignment (blocked version).");
	
	// types
	using SR = typename std::conditional<
		std::is_same<T_OUT, CommonKmerLight>::value,
		KmerIntersectLight<MatrixEntry, T_OUT>,
		KmerIntersectLoc<MatrixEntry, T_OUT>
		>::type;
	using T_SpMat = typename PSpMat<T_OUT>::MPI_DCCols;
	using DER_IN = typename combblas::SpDCCols<uint64_t, MatrixEntry>;

	parops->tp->start_timer("seq-comm-wait");
	
	// wait till receive block seqs
	if (!dfd->is_ready())
	  	dfd->bl_wait();

	parops->tp->stop_timer("seq-comm-wait");

	
	parops->tp->start_timer("sparse");
	parops->tp->start_timer("block-split");

	// block spgemm instance (splits matrices)
	combblas::BlockSpGEMM<uint64_t, MatrixEntry, DER_IN,
						  MatrixEntry, DER_IN>
		bspgemm(A, AT, params.br, params.bc);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("block-split");

	// block matrix and combblas matrix offsets
	vector<pair<uint64_t, uint64_t>> cb_roffsets =
		cbmat_offsets(A.getnrow(), params.br, false);
	vector<pair<uint64_t, uint64_t>> cb_coffsets =
		cbmat_offsets(AT.getncol(), params.bc, true);


	// final similarity matrix els
	vector<uint64_t> rids;
	vector<uint64_t> cids;
	vector<float> vals;
	
	// aligner
	PWAlign *pwa = nullptr;
	if (params.pw_aln != params_t::PwAln::ALN_NONE)
	{
		if (params.pw_aln == params_t::PwAln::ALN_SEQAN_FULL)
			pwa = new SeqanFullAligner(params.gap_open, params.gap_ext);
		else if (params.pw_aln == params_t::PwAln::ALN_SEQAN_XDROP)
			pwa = new SeqanXdropAligner
				(params.gap_open, params.gap_ext, params.klength,
				 params.aln_seqan_xdrop, params.seed_count);
#ifdef ADEPT_ENABLED
	    	else if (params.pw_aln == params_t::PwAln::ALN_ADEPT_GPUBSW)
			pwa = new AdeptBSWAligner(params.gap_open, params.gap_ext);
#endif
		else if (params.pw_aln == params_t::PwAln::ALN_IPUMA)
			pwa = new IPumaAligner(params.gap_open, params.gap_ext, params.klength, params.seed_count);
				
		parops->tp->start_timer("construct-seqs");
		
		// raw fasta data to aligner-specific representation
		pwa->construct_seqs_bl(dfd);

		parops->tp->stop_timer("construct-seqs"); 
	}


	// setup thread data for overlapping
	MultData<MatrixEntry, T_OUT> md(&bspgemm, &params);

	parops->tp->start_timer("mult-align");

	// blocked multiplication and alignment
	uint64_t tot_cnnz = 0;
	bool last_block_aligned = false;
	while (1)
	{
		int bri = md.bri_, bci = md.bci_;
		if (bci == -1)
			bci = 0;		
		parops->info("Block " + to_string(bri) + " " +
						 to_string(bci) + " begin");

		auto t_block = std::chrono::system_clock::now();

		// first but the other blocks will be overlapped with alignments
		if (md.bri_ == 0 && md.bci_ == -1)
			md.computeNext_trg();

		parops->info("Matrix mult complete");

		parops->tp->start_timer("sparse");
		parops->tp->start_timer("prune");		

		auto &C = md.C_[md.curC_];

		uint64_t cnnz = C.getnnz();
		parops->tp->mat_stats["c_nnz"] += cnnz;

		uint64_t roffset = md.roffset_, coffset = md.coffset_;
		if (md.partial_block_)	// get rid of lower triangular matrix elems
		{
			auto f_strict_upper =
				[&roffset, &coffset]
				(const tuple<uint64_t, uint64_t, T_OUT> &t)
				{
					return std::get<0>(t)+roffset >= std::get<1>(t)+coffset;
				};
			C.PruneI(f_strict_upper);
		}

		parops->tp->stop_timer("prune");
		parops->tp->stop_timer("sparse");	

		uint64_t cnnz_partial = C.getnnz();
		parops->tp->mat_stats["prune_sym"] += cnnz_partial;

		if (pwa)
			pw_aln_batch_ovlp(dfd, C, *pwa, params, md, bri, bci);

		uint64_t cnnz_aln = C.getnnz();
		parops->info
			("Block " + to_string(bri) + " " + to_string(bci) + 
			 " computed C nnz " + to_string(cnnz) +
			 " before aln C (partial block) " + to_string(cnnz_partial) +
			 " after aln " + to_string(C.getnnz()) +
			 " elim " + 
			 to_string((static_cast<double>(cnnz-cnnz_aln)/
						cnnz)*100.0) + "%");
		tot_cnnz += cnnz_aln;
		parops->tp->mat_stats["sim_thr"] += cnnz_aln;

		if (C.seqptr()->getnnz() != 0)
		{
			parops->tp->start_timer("mat-formation");

			auto dcsc = C.seqptr()->GetDCSC();
			for (uint64_t i = 0; i < dcsc->nzc; ++i)
			{
				for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
				{
					rids.push_back(dcsc->ir[j]+cb_roffsets[bri].first);
					cids.push_back(dcsc->jc[i]+cb_coffsets[bci].first);
					vals.push_back(dcsc->numx[j].score_aln);
				}
			}

			parops->tp->stop_timer("mat-formation");

			#if PASTIS_DBG_LVL > 0
			parops->bytes_alloc += rids.capacity() * sizeof(uint64_t);
			parops->bytes_alloc += cids.capacity() * sizeof(uint64_t);
			parops->bytes_alloc += vals.capacity() * sizeof(float);
			parops->logger->log("approximate memory in usage " +
								gb_str(parops->bytes_alloc));
			#endif
		}

		double t_block_cur =
			ms_t(std::chrono::system_clock::now()-t_block).count();
		parops->info
			("Block " + to_string(bri) + " " + to_string(bci) +
			 " took " + to_string(t_block_cur/1e3) + " seconds");

		if (last_block_aligned || params.br*params.bc==1)
			break;
		if (md.bri_ == params.br-1 && md.bci_ == params.bc-1)
			last_block_aligned = true;
	}
	
	parops->tp->stop_timer("mult-align");

	// parops->tp->stop_timer("total");
	// string stats = parops->tp->to_string();
	// if (parops->g_rank == 0)
	// 	std::cout << stats << endl;
	// parops->tp->start_timer("total");
	
	parops->tp->start_timer("sparse");
	parops->tp->start_timer("sim-mat-assemble");

	combblas::FullyDistVec<uint64_t, uint64_t> drows(rids, parops->grid);
	combblas::FullyDistVec<uint64_t, uint64_t> dcols(cids, parops->grid);
	combblas::FullyDistVec<uint64_t, float> dvals(vals, parops->grid);
	combblas::SpParMat<uint64_t, float,
					   combblas::SpDCCols<uint64_t, float>>
		S(A.getnrow(), AT.getncol(), drows, dcols, dvals);

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("sim-mat-assemble");
	parops->tp->start_timer("sim-mat-io");
	
	if (!params.align_file.empty())
		S.ParallelWriteMM(params.align_file, true);

	parops->tp->stop_timer("sim-mat-io");
	
	parops->info("tot nnz in similarity mat " + to_string(S.getnnz()));

	if (pwa)
		delete pwa;
}

	


template <typename T_OUT>
void
mult_aln
(
 	params_t						&params,
   	shared_ptr<DistFastaData>		 dfd,
 	PSpMat<MatrixEntry>::MPI_DCCols &A,
	PSpMat<MatrixEntry>::MPI_DCCols &AT
)
{
	parops->logger->log("Computing SpGEMM and alignment (monolith version).");
	
	// types
	using SR = typename std::conditional<
		std::is_same<T_OUT, CommonKmerLight>::value,
		KmerIntersectLight<MatrixEntry, T_OUT>,
		KmerIntersectLoc<MatrixEntry, T_OUT>
		>::type;
	using T_SpMat = typename PSpMat<T_OUT>::MPI_DCCols;


	// output matrix
	parops->logger->log("Computing SpGEMM AA^T.");
	parops->tp->start_timer("main:(AS)AT");
	T_SpMat C = combblas::Mult_AnXBn_DoubleBuff
		<SR, T_OUT, combblas::SpDCCols<uint64_t, T_OUT>>(A, AT);
	parops->tp->stop_timer("main:(AS)AT");
	C.PrintInfo();


	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc += C.getlocalnnz() * sizeof(uint64_t);
	parops->bytes_alloc += C.getlocalnnz() * sizeof(T_OUT);
	parops->bytes_alloc += C.getlocalcols() * sizeof(uint64_t) * 2;
	parops->logger->log("approximate memory in usage " +
						gb_str(parops->bytes_alloc));
	#endif


	// @OGUZ-TODO symmetricize output matrix when subs are used

	parops->info("C (overlap matrix) load imbalance: " +
				 to_string(C.LoadImbalance()));


	// wait till receive grid seqs
	parops->tp->start_timer("main:dfd->wait");
	if (!dfd->is_ready())
	  	dfd->wait();
	parops->tp->stop_timer("main:dfd->wait");
	

	parops->tp->start_timer("main:sim->align");

	if (params.pw_aln != params_t::PwAln::ALN_NONE)
	{
		PWAlign *pwa = nullptr;
		if (params.pw_aln == params_t::PwAln::ALN_SEQAN_FULL)
			pwa = new SeqanFullAligner(params.gap_open, params.gap_ext);
		else if (params.pw_aln == params_t::PwAln::ALN_SEQAN_XDROP)
			pwa = new SeqanXdropAligner
				(params.gap_open, params.gap_ext, params.klength,
				 params.aln_seqan_xdrop, params.seed_count);
#ifdef ADEPT_ENABLED	
		else if (params.pw_aln == params_t::PwAln::ALN_ADEPT_GPUBSW)
			pwa = new AdeptBSWAligner(params.gap_open, params.gap_ext);
#endif
		else if (params.pw_aln == params_t::PwAln::ALN_IPUMA)
			pwa = new IPumaAligner(params.gap_open, params.gap_ext, params.klength, params.seed_count);

		parops->tp->start_timer("sim:construct_seqs");

		// raw fasta data to aligner-specific representation
		pwa->construct_seqs(dfd);

		parops->tp->stop_timer("sim:construct_seqs");
		
		pw_aln_batch(dfd, C, *pwa, params);

		delete pwa;
	}

	parops->tp->stop_timer("main:sim->align");


	combblas::SpParMat<uint64_t, float,
					   combblas::SpDCCols<uint64_t, float>> S(C);
	if (!params.align_file.empty())
		S.ParallelWriteMM(params.align_file, true);

	parops->info("tot nnz in similarity mat " + to_string(S.getnnz()));
	

	delete params.alph;
}




template <typename T_OUT>
void
compute_sim_mat
(
    params_t &params
)
{
	parops->logger->log("Reading and distributing sequences.");
	string s_tmp;
	
	// Read and distribute sequences
	parops->tp->start_timer("fasta");
	
	shared_ptr<DistFastaData> dfd =
		std::make_shared<DistFastaData>
		(params.input_file.c_str(), params.idx_map_file,
		 params.input_overlap, params.klength, params);

	parops->tp->stop_timer("fasta");

	if (dfd->global_count() != params.seq_count)
	{
		uint64_t final_seq_count = dfd->global_count();
		s_tmp = "\nINFO: Modified sequence count\n";
		s_tmp.append("  Final sequence count: ")
			.append(std::to_string(final_seq_count))
			.append(" (")
			.append(to_string((((params.seq_count-final_seq_count)
								* 100.0) / params.seq_count)))
			.append("% removed)");
		params.seq_count = dfd->global_count();
		parops->info(s_tmp);
	}


	// generate seq-by-kmer matrix
	parops->tp->start_timer("seq-kmer-mat");
	
	unordered_set<Kmer, Kmer> local_kmers;	
	PSpMat<MatrixEntry>::MPI_DCCols A =
		generate_A(params.seq_count, dfd, params.klength,
				   params.kstride, *params.alph, local_kmers);
	
	parops->info("A (seq-by-kmer matrix) load imbalance: " +
		to_string(A.LoadImbalance()));
	A.PrintInfo();

	parops->tp->stop_timer("seq-kmer-mat");

	parops->tp->start_timer("sparse");
	parops->tp->start_timer("tr");

	// AT
	parops->logger->log("Computing transpose of A.");	
	auto AT = A;
	AT.Transpose();	
	AT.PrintInfo();

	parops->tp->stop_timer("sparse");
	parops->tp->stop_timer("tr");


	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc += AT.getlocalnnz() * sizeof(uint64_t);
	parops->bytes_alloc += AT.getlocalnnz() * sizeof(MatrixEntry);
	parops->bytes_alloc += AT.getlocalcols() * sizeof(uint64_t) * 2;
	#endif

	parops->tp->start_timer("sim-search");

	// multiplication and alignment
	if (params.br == -1)
		mult_aln<T_OUT>(params, dfd, A, AT);
	else
	{
		if (params.pb && params.lb == params_t::LoadBal::LB_TRG)
			mult_aln_bl_trg_ovlp<T_OUT>(params, dfd, A, AT);
		else if (params.pb && params.lb == params_t::LoadBal::LB_IDX)
			mult_aln_bl_idx_ovlp<T_OUT>(params, dfd, A, AT);
		else if (!params.pb && params.lb == params_t::LoadBal::LB_TRG)
			mult_aln_bl_trg<T_OUT>(params, dfd, A, AT);
		else if (!params.pb && params.lb == params_t::LoadBal::LB_IDX)
			mult_aln_bl_idx<T_OUT>(params, dfd, A, AT);
	}

	parops->tp->stop_timer("sim-search");
}



template
void
compute_sim_mat<CommonKmerLight> (params_t &);

template
void
compute_sim_mat<CommonKmerLoc> (params_t &);

}

