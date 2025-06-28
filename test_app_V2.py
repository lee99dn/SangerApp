import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Any
import io
import warnings
from datetime import datetime
import time
from dataclasses import dataclass

# BioPython imports with error handling
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Align import PairwiseAligner
    from Bio.Blast import NCBIWWW, NCBIXML

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


@dataclass
class SequenceData:
    """Data class to hold sequence information."""
    sequence: Seq
    quality_scores: List[int]
    annotations: Dict[str, Any]
    name: str


@dataclass
class TrimResult:
    """Data class to hold trimming results."""
    trimmed_seq: Seq
    trimmed_qual: List[int]
    start_pos: int
    end_pos: int


@dataclass
class AlignmentResult:
    """Data class to hold alignment results."""
    consensus: Seq
    score: float
    alignment_string: str
    stats: Dict[str, Any]


@dataclass
class BlastHit:
    """Data class to hold BLAST hit information."""
    accession: str
    description: str
    score: float
    e_value: float
    identity: float
    coverage: float
    alignment_length: int


class SequenceQualityTrimmer:
    """Handles sequence quality trimming operations."""

    @staticmethod
    def trim_by_quality(seq: Seq, qual: List[int], cutoff: int = 20,
                        window_size: int = 10,
                        score_margin: float = 1.5,
                        warn_callback: Optional[callable] = None) -> TrimResult:
        """
        Hybrid Sanger trimming logic:
        - Finds the first acceptable window (sliding window style).
        - Also finds the best-scoring window above the cutoff.
        - If best is significantly better than first (by score_margin), use it. Else, use the first one.
        """

        # Input validation
        if len(seq) != len(qual):
            raise ValueError(f"Sequence length ({len(seq)}) doesn't match quality length ({len(qual)})")

        if not qual or len(qual) == 0:
            return TrimResult(seq, qual, 0, len(seq) - 1)

        # Handle case where window is larger than sequence
        if window_size > len(qual):
            avg_qual = np.mean(qual)
            if avg_qual >= cutoff:
                return TrimResult(seq, qual, 0, len(seq) - 1)
            else:
                msg = f"Sequence too short for window size {window_size}. Average quality {avg_qual:.1f} below cutoff {cutoff}."
                if warn_callback:
                    warn_callback(msg)
                else:
                    warnings.warn(msg)
                return TrimResult(seq, qual, 0, len(seq) - 1)

        if len(seq) != len(qual):
            raise ValueError("Sequence and quality score lengths do not match.")

    # Handle edge cases
        if len(seq) < window_size:
            return TrimResult(seq, qual, 0, len(seq) - 1)

        first_start = None
        first_mean = None
        best_start = None
        best_mean = -np.inf

        for i in range(len(qual) - window_size + 1):
            window = qual[i:i + window_size]
            window_mean = np.mean(window)

            if window_mean >= cutoff:
                if first_start is None:
                    first_start = i
                    first_mean = window_mean
                if window_mean > best_mean:
                    best_mean = window_mean
                    best_start = i

        # Decide which to use - FIXED LOGIC
        if first_start is None:
        # No window passed the cutoff, return full sequence
            return TrimResult(seq, qual, 0, len(seq) - 1)
    
    # Fixed logic: When score_margin is 0, ALWAYS use first_start
    # When score_margin > 0, only use best if it's better by MORE than the margin
        if score_margin <= 0:
            start = first_start
        else:
            # Only use best if it's significantly better than first by MORE than the margin
            use_best = (best_start is not None and 
                       first_mean is not None and 
                       (best_mean - first_mean) > score_margin)
            start = best_start if use_best else first_start
        
        trimmed_seq = seq[start:]
        trimmed_qual = qual[start:]

        return TrimResult(trimmed_seq, trimmed_qual, start, len(seq) - 1)


class SequenceAligner:
    """Handles sequence alignment operations."""

    @staticmethod
    def _extract_alignment_strings(alignment, seq1: str, seq2: str) -> Tuple[str, str]:
        """Extract aligned strings from BioPython alignment object."""
        aligned_seqA = []
        aligned_seqB = []

        last1 = last2 = 0
        for (start1, end1), (start2, end2) in zip(alignment.aligned[0], alignment.aligned[1]):
            # Add gaps for unaligned regions
            if start1 > last1:
                aligned_seqA.append(seq1[last1:start1])
                aligned_seqB.append('-' * (start1 - last1))
            if start2 > last2:
                aligned_seqA.append('-' * (start2 - last2))
                aligned_seqB.append(seq2[last2:start2])

            # Add aligned region
            aligned_seqA.append(seq1[start1:end1])
            aligned_seqB.append(seq2[start2:end2])

            last1 = end1
            last2 = end2

        return ''.join(aligned_seqA), ''.join(aligned_seqB)

    @staticmethod
    def align_sequences(seq1: Seq, seq2: Seq) -> AlignmentResult:
        """Align two sequences and generate consensus."""
        try:
        # EMBOSS-style scoring parameters
            aligner = PairwiseAligner()
            aligner.match_score = 5
            aligner.mismatch_score = -4
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -0.5
            aligner.mode = 'global'
        # Perform alignment
            alignments = aligner.align(seq1, seq2)
            best_alignment = alignments[0]
        # Extract aligned sequences
            aligned_seqA, aligned_seqB = SequenceAligner._extract_alignment_strings(
                best_alignment, str(seq1), str(seq2)
            )
        # Generate consensus
            consensus_chars, matches, mismatches, gaps = SequenceAligner._generate_sanger_consensus(
                    aligned_seqA, aligned_seqB)

            consensus = Seq(''.join(consensus_chars))
        # Calculate scores and statistics
            max_possible_score = max(len(seq1), len(seq2)) * aligner.match_score
            normalized_score = best_alignment.score / max_possible_score if max_possible_score > 0 else 0
            total_positions = len(aligned_seqA)
            stats = {
                'matches': matches,
                'mismatches': mismatches,
                'gaps': gaps,
                'identity': matches / total_positions * 100 if total_positions > 0 else 0,
                'coverage': (total_positions - gaps) / total_positions * 100 if total_positions > 0 else 0,
                'raw_score': best_alignment.score,
                'max_possible_score': max_possible_score
            }
            return AlignmentResult(consensus, normalized_score, str(best_alignment), stats)
        except Exception as e:
            st.error(f"Error during alignment: {e}")
            return AlignmentResult(Seq(""), 0.0, "", {})


    @staticmethod
    def _generate_sanger_consensus(aligned_seqA: str, aligned_seqB: str):
        """Generate consensus using Sanger-specific logic (trust reverse after first match)."""
        consensus_chars = []
        matches = mismatches = gaps = 0
        first_match_found = False

        for a, b in zip(aligned_seqA, aligned_seqB):
            if a == b and a != '-':
                consensus_chars.append(a)
                matches += 1
                if not first_match_found:
                    first_match_found = True
            elif a == '-':
                consensus_chars.append(b)
                gaps += 1
            elif b == '-':
                consensus_chars.append(a)
                gaps += 1
            else:
            # Handle mismatches - after first match, trust reverse read (aligned_seqB)
                if first_match_found:
                    consensus_chars.append(b)  # Trust reverse read after first match
                else:
                    consensus_chars.append(a)  # Use forward read before first match
                mismatches += 1

        return consensus_chars, matches, mismatches, gaps

        
class BlastAnalyzer:
    """Handles BLAST operations."""

    @staticmethod
    def run_blast_search(sequence: Seq, database: str = "nt", program: str = "blastn",
                         max_hits: int = 10) -> List[BlastHit]:
        """Run BLAST search against NCBI database."""
        try:
            st.info("🔍 Running BLAST search... This may take a few minutes.")
            progress_bar = st.progress(0)
            status_text = st.empty()

            # Submit BLAST job
            status_text.text("Submitting BLAST job...")
            progress_bar.progress(20)

            # BLAST request
            result_handle = NCBIWWW.qblast(
                program=program,
                database=database,
                sequence=str(sequence),
                hitlist_size=max_hits,
                expect=10,
                word_size=11 if program == "blastn" else 3
            )

            # Parse and limit hits
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            hits = []
            for alignment in blast_record.alignments:
                best_hsp = max(alignment.hsps, key=lambda h: h.score)

                hit = BlastHit(
                    accession=alignment.accession,
                    description=alignment.title[:100] + "..." if len(alignment.title) > 100 else alignment.title,
                    score=best_hsp.score,
                    e_value=best_hsp.expect,
                    identity=(best_hsp.identities / best_hsp.align_length) * 100,
                    coverage=(best_hsp.align_length / len(sequence)) * 100,
                    alignment_length=best_hsp.align_length
                )
                hits.append(hit)

                if len(hits) >= max_hits:
                    break

            progress_bar.progress(100)
            status_text.text("BLAST search completed!")
            time.sleep(1)
            progress_bar.empty()
            status_text.empty()

            return hits

        except Exception as e:
            st.error(f"BLAST search failed: {e}")
            return []

    @staticmethod
    def display_blast_results(hits: List[BlastHit]) -> pd.DataFrame:
        """Display BLAST results in a formatted table and return DataFrame for download."""
        if not hits:
            st.warning("No BLAST hits found.")
            return pd.DataFrame()

        # Convert hits to DataFrame
        hit_data = []
        for i, hit in enumerate(hits, 1):
            hit_data.append({
                'Rank': i,
                'Accession': hit.accession,
                'Description': hit.description,
                'Score': f"{hit.score:.1f}",
                'E-value': f"{hit.e_value:.2e}",
                'Identity (%)': f"{hit.identity:.1f}",
                'Coverage (%)': f"{hit.coverage:.1f}",
                'Alignment Length': hit.alignment_length
            })

        df = pd.DataFrame(hit_data)
        top_50_df = df.head(50)
        st.info("TOP 50")
        st.dataframe(top_50_df, use_container_width=True, hide_index=True)

        # Summary statistics
        if hits:
            best_hit = hits[0]
            st.info(
                f"🎯 Best match: {best_hit.description} (Identity: {best_hit.identity:.1f}%, E-value: {best_hit.e_value:.2e})")

        return df


class FileHandler:
    """Handles file I/O operations."""

    @staticmethod
    def read_ab1_file(file_content: bytes) -> Optional[SequenceData]:
        """Read AB1 file and extract sequence data."""
        try:
            file_obj = io.BytesIO(file_content)
            record = SeqIO.read(file_obj, "abi")

            return SequenceData(
                sequence=record.seq,
                quality_scores=record.letter_annotations.get("phred_quality", []),
                annotations=record.annotations,
                name=record.id or "Unknown"
            )
        except Exception as e:
            st.error(f"Error reading AB1 file: {e}")
            return None


class SequenceAnalyzer:
    """Main sequence analysis utilities."""

    @staticmethod
    def calculate_stats(seq: Seq, qual: Optional[List[int]] = None) -> Dict[str, Any]:
        """Calculate sequence statistics."""
        seq_str = str(seq)
        stats = {
            'Length': len(seq_str),
            'A': seq_str.count('A'),
            'T': seq_str.count('T'),
            'G': seq_str.count('G'),
            'C': seq_str.count('C'),
            'N': seq_str.count('N'),
            'GC Content (%)': round((seq_str.count('G') + seq_str.count('C')) / len(seq_str) * 100, 2) if len(
                seq_str) > 0 else 0
        }

        if qual:
            stats.update({
                'Min Quality': min(qual),
                'Max Quality': max(qual),
                'Mean Quality': round(np.mean(qual), 2),
                'Median Quality': round(np.median(qual), 2)
            })

        return stats


class Visualizer:
    """Handles visualization operations."""

    @staticmethod
    def plot_quality_scores(quality_scores: List[int], title: str = "Quality Scores",
                            trimmed_region: Optional[Tuple[int, int]] = None):
        """Plot quality scores with optional trimmed region highlighting."""
        fig, ax = plt.subplots(figsize=(12, 6))

        positions = range(len(quality_scores))
        ax.plot(positions, quality_scores, 'b-', alpha=0.7, linewidth=1)
        ax.fill_between(positions, quality_scores, alpha=0.3)

        # Highlight trimmed region if provided
        if trimmed_region:
            start, end = trimmed_region
            ax.axvspan(start, end, alpha=0.2, color='green', label=f'Kept region ({start}-{end})')
            if start > 0:
                ax.axvspan(0, start, alpha=0.2, color='red', label=f'Trimmed start (0-{start - 1})')
            if end < len(quality_scores) - 1:
                ax.axvspan(end + 1, len(quality_scores) - 1, alpha=0.2, color='red',
                           label=f'Trimmed end ({end + 1}-{len(quality_scores) - 1})')

        ax.set_xlabel('Position')
        ax.set_ylabel('Quality Score')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        if trimmed_region:
            ax.legend()

        return fig


class SangerAnalysisApp:
    """Main application class."""

    def __init__(self):
        self.trimmer = SequenceQualityTrimmer()
        self.aligner = SequenceAligner()
        self.blast_analyzer = BlastAnalyzer()
        self.file_handler = FileHandler()
        self.analyzer = SequenceAnalyzer()
        self.visualizer = Visualizer()

    def setup_page(self):
        """Setup Streamlit page configuration."""
        st.set_page_config(
            page_title="Enhanced Sanger Sequencing Analysis",
            page_icon="🧬",
            layout="wide"
        )

        st.title("🧬 Enhanced Sanger Sequencing Analysis")
        st.markdown(
            "Upload AB1 chromatogram files for forward and reverse reads to generate consensus sequences and perform BLAST analysis")

        if not BIOPYTHON_AVAILABLE:
            st.error("BioPython is required for this application. Please install it with: pip install biopython")
            st.stop()

    def setup_sidebar(self) -> Dict[str, Any]:
        """Setup sidebar parameters."""
        st.sidebar.header("Analysis Parameters")

        # Quality trimming parameters
        st.sidebar.subheader("Quality Trimming")
        quality_cutoff = st.sidebar.slider("Quality Cutoff", 10, 40, 20, help="Minimum quality score threshold")
        window_size = st.sidebar.slider("Window Size", 2, 30, 10, help="Size of sliding window for quality assessment")
        score_margin = st.sidebar.slider("Score Margin", 0.0, 5.0, 1.5, step=0.1, help="How much better the best window must be to override the first acceptable window")


        # BLAST parameters
        st.sidebar.subheader("BLAST Analysis")
        run_blast = st.sidebar.checkbox("Run BLAST Analysis", False, help="Perform BLAST search on consensus sequence")
        blast_database = st.sidebar.selectbox("BLAST Database", ["nt", "nr", "16S ribosomal RNA"], help="NCBI database to search")
        blast_program = st.sidebar.selectbox("BLAST Program", ["blastn", "blastx"], help="BLAST program to use")
        max_hits = st.sidebar.slider("Max BLAST Hits", 5, 500, 10, help="Maximum number of BLAST hits to retrieve")

        return {
            'quality_cutoff': quality_cutoff,
            'window_size': window_size,
            'score_margin': score_margin,
            'run_blast': run_blast,
            'blast_database': blast_database,
            'blast_program': blast_program,
            'max_hits': max_hits
        }

    def handle_file_upload(self) -> Tuple[Optional[SequenceData], Optional[SequenceData]]:
        """Handle file upload and reading."""
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Forward Read (AB1)")
            forward_file = st.file_uploader(
                "Upload Forward AB1 file",
                type=['ab1', 'abi'],
                key="forward"
            )

        with col2:
            st.subheader("Reverse Read (AB1)")
            reverse_file = st.file_uploader(
                "Upload Reverse AB1 file",
                type=['ab1', 'abi'],
                key="reverse"
            )

        forward_data = reverse_data = None

        if forward_file is not None and reverse_file is not None:
            try:
                forward_data = self.file_handler.read_ab1_file(forward_file.read())
                reverse_data = self.file_handler.read_ab1_file(reverse_file.read())

                if forward_data is None or reverse_data is None:
                    st.error("Failed to read AB1 files")
                    return None, None

                st.success("✅ AB1 files loaded successfully!")

            except Exception as e:
                st.error(f"An error occurred during file reading: {e}")
                return None, None

        return forward_data, reverse_data

    def display_raw_sequence_info(self, forward_data: SequenceData, reverse_data: SequenceData):
        """Display raw sequence information."""
        st.header("📊 Raw Sequence Information")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Forward Read")
            forward_stats = self.analyzer.calculate_stats(forward_data.sequence, forward_data.quality_scores)
            st.dataframe(pd.DataFrame([forward_stats]).T)

            fig_forward = self.visualizer.plot_quality_scores(forward_data.quality_scores,
                                                              "Forward Read Quality Scores")
            st.pyplot(fig_forward)
            plt.close()

        with col2:
            st.subheader("Reverse Read (Original)")
            reverse_stats = self.analyzer.calculate_stats(reverse_data.sequence, reverse_data.quality_scores)
            st.dataframe(pd.DataFrame([reverse_stats]).T)

            fig_reverse = self.visualizer.plot_quality_scores(reverse_data.quality_scores,
                                                              "Reverse Read Quality Scores (Original)")
            st.pyplot(fig_reverse)
            plt.close()

    def perform_quality_trimming(self, forward_data: SequenceData, reverse_data: SequenceData,
                                 params: Dict[str, Any]) -> Tuple[TrimResult, TrimResult, Seq]:
        """Perform quality trimming on sequences."""
        st.header("✂️ Quality Trimming")

        # Trim sequences
        forward_trimmed = self.trimmer.trim_by_quality(
            forward_data.sequence, forward_data.quality_scores,
            cutoff=params['quality_cutoff'], window_size=params['window_size'],
            score_margin=params['score_margin'], warn_callback=st.warning
        )

        reverse_trimmed = self.trimmer.trim_by_quality(
            reverse_data.sequence, reverse_data.quality_scores,
            cutoff=params['quality_cutoff'], window_size=params['window_size'],
           score_margin=params['score_margin'], warn_callback=st.warning
        )

        # Apply reverse complement to reverse sequence
        reverse_trimmed_rc = reverse_trimmed.trimmed_seq.reverse_complement()

        # Display results
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Forward Read - Trimmed")
            forward_stats = self.analyzer.calculate_stats(forward_trimmed.trimmed_seq, forward_trimmed.trimmed_qual)
            st.dataframe(pd.DataFrame([forward_stats]).T)

            fig_forward_trim = self.visualizer.plot_quality_scores(
                forward_data.quality_scores,
                "Forward Read - Quality Trimming",
                trimmed_region=(forward_trimmed.start_pos, forward_trimmed.end_pos)
            )
            st.pyplot(fig_forward_trim)
            plt.close()

        with col2:
            st.subheader("Reverse Read - Trimmed (Original orientation)")
            reverse_stats = self.analyzer.calculate_stats(reverse_trimmed.trimmed_seq, reverse_trimmed.trimmed_qual)
            st.dataframe(pd.DataFrame([reverse_stats]).T)

            fig_reverse_trim = self.visualizer.plot_quality_scores(
                reverse_data.quality_scores,
                "Reverse Read - Quality Trimming (Original)",
                trimmed_region=(reverse_trimmed.start_pos, reverse_trimmed.end_pos)
            )
            st.pyplot(fig_reverse_trim)
            plt.close()

        # Show reverse complement
        st.subheader("🔄 Reverse Read After Reverse Complement")
        st.info(
            "This is the reverse-complemented and trimmed reverse sequence that will be used for consensus generation.")

        col1, col2 = st.columns(2)
        with col1:
            st.text("Reverse Trimmed (Original):")
            st.code(str(reverse_trimmed.trimmed_seq)[:100] + ("..." if len(reverse_trimmed.trimmed_seq) > 100 else ""))
        with col2:
            st.text("Reverse Trimmed (Reverse Complement):")
            st.code(str(reverse_trimmed_rc)[:100] + ("..." if len(reverse_trimmed_rc) > 100 else ""))

        return forward_trimmed, reverse_trimmed, reverse_trimmed_rc

    def generate_consensus(self, forward_trimmed: TrimResult, reverse_trimmed_rc: Seq,
                           params: Dict[str, Any]) -> Optional[AlignmentResult]:
        """Generate consensus sequence."""
        st.header("🔗 Consensus Generation")

        if len(forward_trimmed.trimmed_seq) > 0 and len(reverse_trimmed_rc) > 0:
            alignment_result = self.aligner.align_sequences(
                forward_trimmed.trimmed_seq, reverse_trimmed_rc
            )

            if len(alignment_result.consensus) > 0:
                # Display results
                col1, col2 = st.columns([2, 1])

                with col1:
                    st.subheader("Consensus Sequence")
                    st.text_area(
                        "Consensus (FASTA format)",
                        f">Consensus_Sequence_{datetime.now().strftime('%Y%m%d_%H%M%S')}\n{alignment_result.consensus}",
                        height=150
                    )

                    consensus_stats = self.analyzer.calculate_stats(alignment_result.consensus)
                    st.subheader("Consensus Statistics")
                    st.dataframe(pd.DataFrame([consensus_stats]).T)

                with col2:
                    st.subheader("Alignment Statistics")
                    alignment_data = {
                        'Metric': ['Matches', 'Mismatches', 'Gaps', 'Identity (%)', 'Coverage (%)', 'Alignment Score'],
                        'Value': [
                            alignment_result.stats['matches'],
                            alignment_result.stats['mismatches'],
                            alignment_result.stats['gaps'],
                            f"{alignment_result.stats['identity']:.2f}%",
                            f"{alignment_result.stats['coverage']:.2f}%",
                            f"{alignment_result.score:.3f}"
                        ]
                    }
                    st.dataframe(pd.DataFrame(alignment_data), hide_index=True)

                return alignment_result
            else:
                st.error("Failed to generate consensus sequence")
                return None
        else:
            st.error("Trimmed sequences are too short for consensus generation")
            return None

    def run_blast_analysis(self, consensus: Seq, params: Dict[str, Any]) -> Optional[pd.DataFrame]:
        """Run BLAST analysis on consensus sequence."""
        if params['run_blast'] and len(consensus) > 0:
            st.header("🎯 BLAST Analysis")

            # Minimum sequence length check
            if len(consensus) < 20:
                st.warning("Consensus sequence too short for meaningful BLAST analysis (minimum 20 bp recommended)")
                return None

            # Create a unique key for this BLAST search based on sequence and parameters
            blast_key = f"blast_{hash(str(consensus))}_{params['blast_database']}_{params['blast_program']}_{params['max_hits']}"

            # Initialize session state for BLAST results if not exists
            if 'blast_results' not in st.session_state:
                st.session_state.blast_results = {}

            # Check if we already have results for this exact search
            if blast_key in st.session_state.blast_results:
                st.info("Using cached BLAST results. Click 'Clear BLAST Cache' below to run a new search.")
                blast_hits = st.session_state.blast_results[blast_key]
            else:
                # Run new BLAST search
                blast_hits = self.blast_analyzer.run_blast_search(
                    consensus,
                    database=params['blast_database'],
                    program=params['blast_program'],
                    max_hits=params['max_hits']
                )

                # Cache the results
                if blast_hits:  # Only cache if we got results
                    st.session_state.blast_results[blast_key] = blast_hits

            if blast_hits:
                st.subheader("BLAST Results")
                blast_df = self.blast_analyzer.display_blast_results(blast_hits)

                # Add a button to clear BLAST cache if user wants to run a fresh search
                if st.button("🔄 Clear BLAST Cache & Run New Search"):
                    if blast_key in st.session_state.blast_results:
                        del st.session_state.blast_results[blast_key]
                    st.rerun()

                return blast_df
            else:
                st.warning("No BLAST hits found or search failed.")
                return None
        return None

    def provide_downloads(self, alignment_result: AlignmentResult, forward_trimmed: TrimResult,
                          reverse_trimmed_rc: Seq, params: Dict[str, Any], blast_df: Optional[pd.DataFrame] = None):
        """Provide download options for results."""
        st.header("💾 Download Results")

        # Determine number of columns based on whether BLAST results are available
        if blast_df is not None and not blast_df.empty:
            col1, col2, col3, col4 = st.columns(4)
        else:
            col1, col2, col3 = st.columns(3)
            col4 = None

        with col1:
            # Consensus sequence download
            consensus_fasta = f">Consensus_Sequence_{datetime.now().strftime('%Y%m%d_%H%M%S')}\n{alignment_result.consensus}"
            st.download_button(
                label="Download Consensus (FASTA)",
                data=consensus_fasta,
                file_name=f"consensus_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
                mime="text/plain"
            )

        with col2:
            # Statistics report
            report_data = {
                'Consensus_Length': len(alignment_result.consensus),
                'Alignment_Identity': f"{alignment_result.stats['identity']:.2f}%",
                'Alignment_Coverage': f"{alignment_result.stats['coverage']:.2f}%",
                'Quality_Cutoff': params['quality_cutoff'],
                'Window_Size': params['window_size'],
            }
            report_df = pd.DataFrame([report_data]).T
            report_csv = report_df.to_csv()

            st.download_button(
                label="Download Statistics (CSV)",
                data=report_csv,
                file_name=f"sanger_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )

        with col3:
            # Trimmed sequences
            trimmed_fasta = f">Forward_Trimmed_{datetime.now().strftime('%Y%m%d_%H%M%S')}\n{forward_trimmed.trimmed_seq}\n>Reverse_Trimmed_RC_{datetime.now().strftime('%Y%m%d_%H%M%S')}\n{reverse_trimmed_rc}"
            st.download_button(
                label="Download Trimmed Sequences",
                data=trimmed_fasta,
                file_name=f"trimmed_sequences_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
                mime="text/plain"
            )

        if col4 is not None and blast_df is not None and not blast_df.empty:
            with col4:
                # BLAST results download
                blast_csv = blast_df.to_csv(index=False)
                st.download_button(
                    label="Download BLAST Results (CSV)",
                    data=blast_csv,
                    file_name=f"blast_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv"
                )

                # Also provide BLAST results in tab-delimited format for further analysis
                blast_tsv = blast_df.to_csv(index=False, sep='\t')
                st.download_button(
                    label="Download BLAST Results (TSV)",
                    data=blast_tsv,
                    file_name=f"blast_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv",
                    mime="text/tab-separated-values"
                )

    def run(self):
        """Main application entry point."""
        self.setup_page()
        params = self.setup_sidebar()

        forward_data, reverse_data = self.handle_file_upload()

        if forward_data is not None and reverse_data is not None:
            try:
                # Display raw sequence information
                self.display_raw_sequence_info(forward_data, reverse_data)

                # Perform quality trimming
                forward_trimmed, reverse_trimmed, reverse_trimmed_rc = self.perform_quality_trimming(
                    forward_data, reverse_data, params
                )

                # Generate consensus
                alignment_result = self.generate_consensus(forward_trimmed, reverse_trimmed_rc, params)

                if alignment_result:
                    # Run BLAST analysis if requested
                    blast_df = self.run_blast_analysis(alignment_result.consensus, params)

                    # Provide download options
                    self.provide_downloads(alignment_result, forward_trimmed, reverse_trimmed_rc, params, blast_df)

            except Exception as e:
                st.error(f"An error occurred during analysis: {e}")
                st.exception(e)
        else:
            st.info("👆 Please upload both forward and reverse AB1 files to begin analysis")


def main():
    """Main function to run the application."""
    app = SangerAnalysisApp()
    app.run()


if __name__ == "__main__":
    main()
