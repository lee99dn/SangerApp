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
import re
import matplotlib.ticker as ticker

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
    limit_exceeded: bool = False


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
    scientific_name: str
    score: float
    e_value: float
    identity: float
    coverage: float
    alignment_length: int
    query_length: int
    subject_length: int


class SequenceQualityTrimmer:
    """Handles sequence quality trimming operations."""

    @staticmethod
    def trim_by_quality(seq: Seq, qual: List[int], cutoff: int = 20,
                        window_size: int = 10,
                        adjustment: int = 0,
                        read_name: str = "Sequence",
                        warn_callback: Optional[callable] = None) -> TrimResult:
        
        limit_hit = False

        # Input validation
        if len(seq) != len(qual):
            raise ValueError(f"Sequence length ({len(seq)}) doesn't match quality length ({len(qual)})")

        if not qual or len(qual) == 0:
            return TrimResult(seq, qual, 0, len(seq) - 1, limit_hit)

        if window_size > len(qual):
            avg_qual = np.mean(qual)
            if avg_qual >= cutoff:
                return TrimResult(seq, qual, 0, len(seq) - 1, limit_hit)
            else:
                msg = f"[{read_name}] Sequence too short for window size {window_size}. Average quality {avg_qual:.1f} below cutoff {cutoff}."
                if warn_callback:
                    warn_callback(msg)
                else:
                    warnings.warn(msg)
                return TrimResult(seq, qual, 0, len(seq) - 1, limit_hit)

        if len(seq) < window_size:
            return TrimResult(seq, qual, 0, len(seq) - 1, limit_hit)

        # --- 1. Find 5' Start (Forward Scan) ---
        first_start = None

        for i in range(len(qual) - window_size + 1):
            window = qual[i:i + window_size]
            window_mean = np.mean(window)
            if window_mean >= cutoff:
                first_start = i
                # Limit 5' trim to 30% independently
                if first_start > 0.3 * len(seq):
                   first_start = int(0.3 * len(seq))
                   limit_hit = True
                   msg = f"[{read_name}] 5' trimming exceeds 30% limit. Adjusted start position to bp {first_start + 1}."
                   if warn_callback:
                       warn_callback(msg)
                break 

        if first_start is None:
            msg = f"[{read_name}] No window of size {window_size} found with quality >= {cutoff}. Returning full sequence."
            if warn_callback:
                warn_callback(msg)
            else:
                warnings.warn(msg)
            return TrimResult(seq, qual, 0, len(seq) - 1, limit_hit)

        # --- 2. Find 3' End (Backward Scan) ---
        end = len(qual) - 1
        # backward scan
        for i in range(len(qual) - window_size, first_start - 1, -1):
            window = qual[i:i + window_size]
            if np.mean(window) >= cutoff:
                #Final window instead
                end = i + window_size - 1
                break
                
        end = max(first_start, end)
        end = min(len(seq) - 1, end)

        # --- 3. Limit 3' end trim to 30% independently ---
        max_trim_per_end = int(0.3 * len(seq))
        current_3_trim = (len(seq) - 1) - end
        
        if current_3_trim > max_trim_per_end:
            # Re-adjust end position to respect the 30% limit for the 3' end
            new_end = (len(seq) - 1) - max_trim_per_end
            if new_end > end:
                end = new_end
                limit_hit = True
                msg = f"[{read_name}] 3' trimming exceeds 30% limit. Adjusted 3' trim position to bp {end + 1}."
                if warn_callback:
                    warn_callback(msg)

        # --- 4. Apply Manual Adjustment (5' end only) ---
        start = first_start + adjustment
        start = max(0, start)
        start = min(start, end)

        trimmed_seq = seq[start:end + 1]
        trimmed_qual = qual[start:end + 1]

        return TrimResult(trimmed_seq, trimmed_qual, start, end, limit_hit)

class SequenceAligner:
    """Handles sequence alignment and consensus generation."""

    @staticmethod
    def reverse_complement(seq: Seq) -> Seq:
        """Return the reverse complement of a DNA sequence."""
        return seq.reverse_complement()

    @staticmethod
    def _extract_alignment_strings(alignment, seq1: str, seq2: str) -> Tuple[str, str]:
        """Extract aligned strings from BioPython alignment object, including trailing ends."""
        aligned_seqA = []
        aligned_seqB = []

        last1 = last2 = 0
        for (start1, end1), (start2, end2) in zip(alignment.aligned[0], alignment.aligned[1]):
            if start1 > last1:
                aligned_seqA.append(seq1[last1:start1])
                aligned_seqB.append('-' * (start1 - last1))
            if start2 > last2:
                aligned_seqA.append('-' * (start2 - last2))
                aligned_seqB.append(seq2[last2:start2])

            aligned_seqA.append(seq1[start1:end1])
            aligned_seqB.append(seq2[start2:end2])

            last1 = end1
            last2 = end2

        tail1 = seq1[last1:]
        tail2 = seq2[last2:]
        if len(tail1) > len(tail2):
            tail2 += '-' * (len(tail1) - len(tail2))
        elif len(tail2) > len(tail1):
            tail1 += '-' * (len(tail2) - len(tail1))

        aligned_seqA.append(tail1)
        aligned_seqB.append(tail2)

        return ''.join(aligned_seqA), ''.join(aligned_seqB)

    @staticmethod
    def _align_quality_to_gapped(aligned_seq: str, raw_qual: List[int]) -> List[int]:
        """Align quality scores to gapped sequences."""
        aligned_qual = []
        q_index = 0
        for base in aligned_seq:
            if base == '-':
                aligned_qual.append(0)
            else:
                aligned_qual.append(raw_qual[q_index])
                q_index += 1
        return aligned_qual

    @staticmethod
    def _generate_consensus_by_quality(aligned_f: str, aligned_r: str,
                                       qual_f: List[int], qual_r: List[int]) -> Tuple[List[str], int, int, int]:
        """Generate consensus trusting phred quality after alignment between forward and reverse."""
        consensus = []
        matches = mismatches = gaps = 0
        length = len(aligned_f)

        for i in range(length):
            a, b = aligned_f[i], aligned_r[i]
            qf, qr = qual_f[i], qual_r[i]

            if a == b:
                consensus.append(a)
                if a != '-':
                    matches += 1
            elif a == '-':
                consensus.append(b)
                gaps += 1
            elif b == '-':
                consensus.append(a)
                gaps += 1
            else:
                # Mismatch case: choose higher quality
                if qf > qr:
                    consensus.append(a)
                elif qr > qf:
                    consensus.append(b)
                # Tie breaking case: quality are equal
                else:
                    consensus.append(a if i < length // 2 else b)
                mismatches += 1
        return consensus, matches, mismatches, gaps

    @staticmethod
    def align_sequences(seq1: Seq, seq2: Seq,
                        qual1: Optional[List[int]] = None,
                        qual2: Optional[List[int]] = None) -> AlignmentResult:
        """Align two sequences with clustalW scoring and generate consensus with quality-based logic."""
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -0.5
        aligner.query_end_gap_score = 0.0
        aligner.target_end_gap_score = 0.0

        alignment = aligner.align(seq1, seq2)[0]
        aligned_f, aligned_r = SequenceAligner._extract_alignment_strings(alignment, str(seq1), str(seq2))

        # Align qualities to gapped sequences
        aligned_qf = SequenceAligner._align_quality_to_gapped(aligned_f, qual1)
        aligned_qr = SequenceAligner._align_quality_to_gapped(aligned_r, qual2)

        # Generate consensus
        consensus_chars, matches, mismatches, gaps = SequenceAligner._generate_consensus_by_quality(
            aligned_f, aligned_r, aligned_qf, aligned_qr
        )
        consensus_seq = Seq(''.join(consensus_chars))

        # Compute stats
        stats = {
            'matches': matches,
            'mismatches': mismatches,
            'gaps': gaps,
        }

        return AlignmentResult(
            consensus=consensus_seq,
            score=alignment.score,
            alignment_string=str(alignment),
            stats=stats
        )


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
                expect=0.000001 if program == "blastn" else 10,
                word_size=11 if program == "blastn" else 3
            )

            # Parse and limit hits
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            hits = []
            for alignment in blast_record.alignments:

                hit = BlastHit(
                    accession=alignment.accession,
                    description=alignment.title[:100] + "..." if len(alignment.title) > 100 else alignment.title,
                    scientific_name=extract_organism_from_blast_title(alignment.title),
                    score=alignment.hsps[0].score,
                    e_value=alignment.hsps[0].expect,
                    identity=round((alignment.hsps[0].identities / alignment.hsps[0].align_length) * 100, 2),
                    coverage=round((alignment.hsps[0].align_length / len(sequence)) * 100, 2) if alignment.hsps[0].align_length < len(sequence) else 100.00,
                    alignment_length=alignment.hsps[0].align_length,
                    query_length=len(str(sequence)),
                    subject_length=int(alignment.length)
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
                'Scientific_name': hit.scientific_name,
                'Score': f"{hit.score:.2f}",
                'E-value': f"{hit.e_value:.2e}",
                'Identity (%)': f"{hit.identity:.2f}",
                'Coverage (%)': f"{hit.coverage:.2f}",
                'Alignment Length': hit.alignment_length,
                'Query Length': hit.query_length,
                'Subject Length': hit.subject_length
            })

        df = pd.DataFrame(hit_data)
        top_50_df = df.head(50)
        st.info("TOP 50")
        st.dataframe(top_50_df, use_container_width=True, hide_index=True)

        # Summary statistics
        if hits:
            best_hit = hits[0]
            st.info(
                f"🎯 Best match: {best_hit.description} (Identity: {best_hit.identity:.2f}%, E-value: {best_hit.e_value:.2e})")

        return df


def extract_organism_from_blast_title(title):
    """
    Extract organism name from various BLAST XML title formats
    """
    title = title.strip()
    if '|' in title:
        parts = title.split('|')
        if len(parts) >= 4:
            description = parts[4].strip() if len(parts) > 4 else parts[-1].strip()
        else:
            description = parts[-1].strip()
    else:
        match = re.match(r'^(\S+)\s+(.+)', title)
        if match:
            description = match.group(2).strip()
        else:
            description = title

    patterns_to_remove = [
        r'\s+16S\s+ribosomal\s+RNA.*$',
        r'\s+16S\s+rRNA.*$',
        r'\s+ribosomal\s+RNA.*$',
        r'\s+rRNA.*$',
        r',?\s+complete\s+sequence.*$',
        r',?\s+partial\s+sequence.*$',
        r',?\s+complete\s+genome.*$',
        r',?\s+chromosome.*$',
        r',?\s+plasmid.*$',
        r'\s+gene.*$'
    ]

    organism_name = description
    for pattern in patterns_to_remove:
        organism_name = re.sub(pattern, '', organism_name, flags=re.IGNORECASE)
        if organism_name != description: 
            break

    return organism_name.strip()


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
            'Metric': ['Length', 'A', 'T', 'G', 'C', 'N', 'GC Content (%)'],
            'Value': [
                len(seq_str),
                seq_str.count('A'),
                seq_str.count('T'),
                seq_str.count('G'),
                seq_str.count('C'),
                seq_str.count('N'),
                round((seq_str.count('G') + seq_str.count('C')) / len(seq_str) * 100, 2) if len(
                    seq_str) > 0 else 0.0
                ]
        }

        if qual:
            stats['Metric'].extend(['Min Quality', 'Max Quality', 'Mean Quality', 'Median Quality'])
            stats['Value'].extend([
                float(min(qual)),
                float(max(qual)),
                round(np.mean(qual), 2),
                round(np.median(qual), 2)
            ])
        return stats


class Visualizer:
    """Handles visualization operations."""

    @staticmethod
    def plot_quality_scores(quality_scores: List[int], title: str = "Quality Scores",
                            trimmed_region: Optional[Tuple[int, int]] = None,
                            window_size: Optional[int] = None, cutoff: Optional[int] = None):
        """Plot quality scores with optional trimmed region highlighting and rolling mean."""
        fig, ax = plt.subplots(figsize=(12, 6))

        # ปรับแกน x ให้เริ่มที่ 1 (1-based indexing)
        positions = range(1, len(quality_scores) + 1)
        
        # Plot Raw Quality
        ax.plot(positions, quality_scores, 'b-', alpha=0.3, linewidth=1, label='Raw Q-Score')

        # Calculate and plot Rolling Mean
        if window_size and cutoff:
            x_mean = []
            y_mean = []
            for i in range(len(quality_scores) - window_size + 1):
                y_mean.append(np.mean(quality_scores[i:i+window_size]))      
                x_mean.append(i + 1 + (window_size - 1) / 2.0)
            ax.plot(x_mean, y_mean, color='darkorange', alpha=1, linewidth=2, 
                    label=f'Window Mean ({window_size}bp)')
            ax.axhline(y=cutoff, color='r', linestyle='--', alpha=0.7, label=f'Cutoff ({cutoff})')

        space = 50
        ax.xaxis.set_major_locator(ticker.MultipleLocator(space))

        # Highlight trimmed region if provided (1-based)
        if trimmed_region:
            start, end = trimmed_region
            start_1b = start + 1
            end_1b = end + 1
            ax.axvspan(start_1b, end_1b, alpha=0.2, color='green', label=f'Kept region ({start_1b}-{end_1b})')
            if start_1b > 1:
                ax.axvspan(1, start_1b - 1, alpha=0.2, color='red', label=f'Trimmed start (1-{start_1b - 1})')
            if end_1b < len(quality_scores):
                ax.axvspan(end_1b + 1, len(quality_scores), alpha=0.2, color='red',
                           label=f'Trimmed end ({end_1b + 1}-{len(quality_scores)})')

        ax.set_xlabel('Base Position (bp)')
        ax.set_ylabel('Quality Score')
        ax.set_title(title)
        
        
        ax.set_xlim(1, len(quality_scores))
        
        
        ax.set_ylim(bottom=0)
        
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')

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
            page_title="SangerApp",
            page_icon="🧬",
            layout="wide"
        )

        st.title("🧬 SangerApp")
        st.markdown(
            "Upload AB1 chromatogram files for forward and reverse reads to generate consensus sequences and perform BLAST analysis")
        url = "https://github.com/lee99dn/SangerAssembly/blob/main/README.md"
        st.markdown("Tutorial >>>> [link](%s)" % url)

        if not BIOPYTHON_AVAILABLE:
            st.error("BioPython is required for this application. Please install it with: pip install biopython")
            st.stop()

    def setup_sidebar(self) -> Dict[str, Any]:
        """Setup sidebar parameters."""
        st.sidebar.header("Analysis Parameters")

        # Quality trimming parameters
        st.sidebar.subheader("Quality Trimming")
        quality_cutoff = st.sidebar.slider("Quality Cutoff", 10, 40, 20, 
                                           help="Minimum quality score threshold based on window mean\n(Total trim cannot exceed 30% of sequence)")
        window_size = st.sidebar.slider("Window Size", 2, 30, 10, help="Size of sliding window for quality assessment")
        adjustment = st.sidebar.slider("Start Adjustment (5' end)", -30, 30, 0, step=1,
                                         help="Adjust the start trim position from the first acceptable window")

        # BLAST parameters
        st.sidebar.subheader("BLAST Analysis")
        run_blast = st.sidebar.checkbox("Run BLAST Analysis", False, help="Perform BLAST search on consensus sequence")
        blast_database = st.sidebar.selectbox("BLAST Database", ["nt", "nr", "Bacteria and Archaea"],
                                              help="NCBI database to search")
        blast_program = st.sidebar.selectbox("BLAST Program", ["blastn", "blastx"], help="BLAST program to use")
        max_hits = st.sidebar.slider("Max BLAST Hits", 5, 500, 10, help="Maximum number of BLAST hits to retrieve")

        # Output parameters
        st.sidebar.subheader("Output Settings")
        use_smart_name = st.sidebar.checkbox("Smart Name", True, 
                                             help="Automatically extract sample name from filename (e.g., gets '1' from ATGC_xxx_1_27F_...)")

        return {
            'quality_cutoff': quality_cutoff,
            'window_size': window_size,
            'adjustment': adjustment,
            'run_blast': run_blast,
            'blast_database': blast_database,
            'blast_program': blast_program,
            'max_hits': max_hits,
            'use_smart_name': use_smart_name
        }

    def handle_file_upload(self) -> Tuple[Optional[SequenceData], Optional[SequenceData], Optional[str]]:
        """Handle file upload and reading, returning filename for smart naming."""
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
        fw_filename = None

        if forward_file is not None and reverse_file is not None:
            try:
                forward_data = self.file_handler.read_ab1_file(forward_file.read())
                reverse_data = self.file_handler.read_ab1_file(reverse_file.read())
                fw_filename = forward_file.name  # เก็บชื่อไฟล์เพื่อนำไปใช้ดึงชื่อตัวอย่าง

                if forward_data is None or reverse_data is None:
                    st.error("Failed to read AB1 files")
                    return None, None, None

                st.success("✅ AB1 files loaded successfully!")

            except Exception as e:
                st.error(f"An error occurred during file reading: {e}")
                return None, None, None

        return forward_data, reverse_data, fw_filename

    def display_raw_sequence_info(self, forward_data: SequenceData, reverse_data: SequenceData):
        """Display raw sequence information."""
        st.header("📊 Raw Sequence Information")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Forward Read")
            st.text("Forward Read (Original):")
            st.code(str(forward_data.sequence)[:100] + ("..." if len(forward_data.sequence) > 100 else ""))
            forward_stats = self.analyzer.calculate_stats(forward_data.sequence, forward_data.quality_scores)
            st.dataframe(pd.DataFrame(forward_stats), hide_index=True)

            fig_forward = self.visualizer.plot_quality_scores(forward_data.quality_scores,
                                                              "Forward Read Quality Scores")
            st.pyplot(fig_forward)
            plt.close()

        with col2:
            st.subheader("Reverse Read (Original)")
            st.text("Reverse Read (Original):")
            st.code(str(reverse_data.sequence)[:100] + ("..." if len(reverse_data.sequence) > 100 else ""))
            reverse_stats = self.analyzer.calculate_stats(reverse_data.sequence, reverse_data.quality_scores)
            st.dataframe(pd.DataFrame(reverse_stats), hide_index=True)

            fig_reverse = self.visualizer.plot_quality_scores(reverse_data.quality_scores,
                                                              "Reverse Read Quality Scores (Original)")
            st.pyplot(fig_reverse)
            plt.close()

    def perform_quality_trimming(self, forward_data: SequenceData, reverse_data: SequenceData,
                                 params: Dict[str, Any]) -> Tuple[TrimResult, TrimResult, Seq, List[int]]:
        """Perform quality trimming on sequences."""
        st.header("✂️ Quality Trimming")

        forward_trimmed = self.trimmer.trim_by_quality(
            forward_data.sequence, forward_data.quality_scores,
            cutoff=params['quality_cutoff'], window_size=params['window_size'],
            adjustment=params['adjustment'], read_name="Forward Read",
            warn_callback=st.warning
        )
        reverse_trimmed = self.trimmer.trim_by_quality(
            reverse_data.sequence, reverse_data.quality_scores,
            cutoff=params['quality_cutoff'], window_size=params['window_size'],
            adjustment=params['adjustment'], read_name="Reverse Read",
            warn_callback=st.warning
        )

        reverse_trimmed_rc = SequenceAligner.reverse_complement(reverse_trimmed.trimmed_seq)
        reverse_qual_rc = reverse_trimmed.trimmed_qual[::-1]
        
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Forward Read - Trimmed")
            forward_stats = self.analyzer.calculate_stats(forward_trimmed.trimmed_seq, forward_trimmed.trimmed_qual)
            st.dataframe(pd.DataFrame(forward_stats), hide_index=True)

            fig_forward_trim = self.visualizer.plot_quality_scores(
                forward_data.quality_scores,
                "Forward Read - Quality Trimming",
                trimmed_region=(forward_trimmed.start_pos, forward_trimmed.end_pos),
                window_size=params['window_size'],
                cutoff=params['quality_cutoff']
            )
            st.pyplot(fig_forward_trim)
            plt.close()
            st.success(f"✅ Kept Region: bp {forward_trimmed.start_pos + 1} to {forward_trimmed.end_pos + 1} "
                       f"(Retained length: {forward_trimmed.end_pos - forward_trimmed.start_pos + 1} bp)")

        with col2:
            st.subheader("Reverse Read - Trimmed (Original orientation)")
            reverse_stats = self.analyzer.calculate_stats(reverse_trimmed.trimmed_seq, reverse_trimmed.trimmed_qual)
            st.dataframe(pd.DataFrame(reverse_stats), hide_index=True)

            fig_reverse_trim = self.visualizer.plot_quality_scores(
                reverse_data.quality_scores,
                "Reverse Read - Quality Trimming (Original)",
                trimmed_region=(reverse_trimmed.start_pos, reverse_trimmed.end_pos),
                window_size=params['window_size'],
                cutoff=params['quality_cutoff']
            )
            st.pyplot(fig_reverse_trim)
            plt.close()
            st.success(f"✅ Kept Region: bp {reverse_trimmed.start_pos + 1} to {reverse_trimmed.end_pos + 1} "
                       f"(Retained length: {reverse_trimmed.end_pos - reverse_trimmed.start_pos + 1} bp)")

        st.subheader("🔄 Reverse Read After Reverse Complement")
        st.info("This is the reverse-complemented and trimmed reverse sequence that will be used for consensus generation.")

        col1, col2 = st.columns(2)
        with col1:
            st.text("Forward Trimmed:")
            st.code(str(forward_trimmed.trimmed_seq)[:100] + ("..." if len(forward_trimmed.trimmed_seq) > 100 else ""))
        with col2:
            st.text("Reverse Trimmed (Original):")
            st.code(str(reverse_trimmed.trimmed_seq)[:100] + ("..." if len(reverse_trimmed.trimmed_seq) > 100 else ""))
            st.text("Reverse Trimmed (Reverse Complement):")
            st.code(str(reverse_trimmed_rc)[:100] + ("..." if len(reverse_trimmed_rc) > 100 else ""))

        return forward_trimmed, reverse_trimmed, reverse_trimmed_rc, reverse_qual_rc

    def generate_consensus(self, forward_trimmed: TrimResult, reverse_trimmed: TrimResult,
                           reverse_trimmed_rc: Seq, reverse_qual_rc: List[int],
                           params: Dict[str, Any], header_id: str) -> Optional[AlignmentResult]:
        """Generate consensus sequence using forward and reverse (rev-comp) trimmed reads."""
        st.header("🔗 Consensus Generation")

        if len(forward_trimmed.trimmed_seq) == 0 or len(reverse_trimmed.trimmed_seq) == 0:
            st.error("Trimmed sequences are too short for consensus generation")
            return None

        alignment_result = self.aligner.align_sequences(
            forward_trimmed.trimmed_seq,
            reverse_trimmed_rc,
            qual1=forward_trimmed.trimmed_qual,
            qual2=reverse_qual_rc
        )

        if len(alignment_result.consensus) == 0:
            st.error("Failed to generate consensus sequence")
            return None

        col1, col2 = st.columns([2, 1])

        with col1:
            st.subheader("Consensus Sequence")
            # ใช้ header_id ที่ส่งเข้ามา
            st.text_area(
                "Consensus (FASTA format)",
                f">{header_id}\n{alignment_result.consensus}",
                height=150
            )

            consensus_stats = self.analyzer.calculate_stats(alignment_result.consensus)
            st.subheader("Consensus Statistics")
            st.dataframe(pd.DataFrame(consensus_stats), hide_index=True)

        with col2:
            st.subheader("Alignment Statistics")
            alignment_data = {
                'Metric': ['Matches', 'Mismatches', 'Gaps', 'Alignment Score'],
                'Value': [
                    alignment_result.stats['matches'],
                    alignment_result.stats['mismatches'],
                    alignment_result.stats['gaps'],
                    f"{alignment_result.score}"
                ]
            }
            st.dataframe(pd.DataFrame(alignment_data), hide_index=True)

        return alignment_result

    def run_blast_analysis(self, consensus: Seq, params: Dict[str, Any]) -> Optional[pd.DataFrame]:
        """Run BLAST analysis on consensus sequence."""
        if params['run_blast'] and len(consensus) > 0:
            st.header("🎯 BLAST Analysis")

            if len(consensus) < 20:
                st.warning("Consensus sequence too short for meaningful BLAST analysis (minimum 20 bp recommended)")
                return None

            blast_key = f"blast_{hash(str(consensus))}_{params['blast_database']}_{params['blast_program']}_{params['max_hits']}"

            if 'blast_results' not in st.session_state:
                st.session_state.blast_results = {}

            if blast_key in st.session_state.blast_results:
                st.info("Using cached BLAST results. Click 'Clear BLAST Cache' below to run a new search.")
                blast_hits = st.session_state.blast_results[blast_key]
            else:
                blast_hits = self.blast_analyzer.run_blast_search(
                    consensus,
                    database=params['blast_database'],
                    program=params['blast_program'],
                    max_hits=params['max_hits']
                )
                if blast_hits:
                    st.session_state.blast_results[blast_key] = blast_hits

            if blast_hits:
                st.subheader("BLAST Results")
                blast_df = self.blast_analyzer.display_blast_results(blast_hits)

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
                          reverse_trimmed: TrimResult, reverse_trimmed_rc: Seq, params: Dict[str, Any], 
                          blast_df: Optional[pd.DataFrame] = None, forward_raw_len: int = 0, reverse_raw_len: int = 0,
                          header_id: str = "Consensus_Sequence"):
        """Provide download options for results."""
        st.header("💾 Download Results")

        if blast_df is not None and not blast_df.empty:
            col1, col2, col3, col4 = st.columns(4)
        else:
            col1, col2, col3 = st.columns(3)
            col4 = None

        with col1:
            # ใช้ header_id ควบคุมข้อมูลและชื่อไฟล์
            consensus_fasta = f">{header_id}\n{alignment_result.consensus}\n"
            st.download_button(
                label="Download Consensus (FASTA)",
                data=consensus_fasta,
                file_name=f"{header_id}.fasta",
                mime="text/plain"
            )

        with col2:
            report_data = {
                "Metric": ["Matches", "Mismatches", "Gaps", "Alignment Score", 'Consensus_Length',
                           'Quality_Cutoff', 'Window_Size', 'Start_Adjustment_5_end', 'Threshold_Type',
                           'Forward_Raw_Length', 'Forward_Trimmed_Length', 'Forward_Kept_Region_bp', 'Forward_Hit_30%_Limit',
                           'Reverse_Raw_Length', 'Reverse_Trimmed_Length', 'Reverse_Kept_Region_bp', 'Reverse_Hit_30%_Limit'],
                "Value": [
                    alignment_result.stats['matches'],
                    alignment_result.stats['mismatches'],
                    alignment_result.stats['gaps'],
                    alignment_result.score,
                    len(alignment_result.consensus),
                    params['quality_cutoff'],
                    params['window_size'],
                    params['adjustment'],
                    "Window Mean",
                    forward_raw_len,
                    len(forward_trimmed.trimmed_seq),
                    f"{forward_trimmed.start_pos + 1} - {forward_trimmed.end_pos + 1}",
                    "Yes" if forward_trimmed.limit_exceeded else "No",
                    reverse_raw_len,
                    len(reverse_trimmed.trimmed_seq),
                    f"{reverse_trimmed.start_pos + 1} - {reverse_trimmed.end_pos + 1}",
                    "Yes" if reverse_trimmed.limit_exceeded else "No"
                ]
            }
            report_df = pd.DataFrame(report_data)
            report_csv = report_df.to_csv(index=False)

            st.download_button(
                label="Download Statistics (CSV)",
                data=report_csv,
                file_name=f"sanger_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )

        with col3:
            trimmed_fasta = f">Forward_Trimmed_{datetime.now().strftime('%Y%m%d_%H%M%S')}\n{forward_trimmed.trimmed_seq}\n>Reverse_Trimmed_RC_{datetime.now().strftime('%Y%m%d_%H%M%S')}\n{reverse_trimmed_rc}"
            st.download_button(
                label="Download Trimmed Sequences",
                data=trimmed_fasta,
                file_name=f"trimmed_sequences_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
                mime="text/plain"
            )

        if col4 is not None and blast_df is not None and not blast_df.empty:
            with col4:
                blast_csv = blast_df.to_csv(index=False)
                st.download_button(
                    label="Download BLAST Results (CSV)",
                    data=blast_csv,
                    file_name=f"blast_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv"
                )

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

        forward_data, reverse_data, fw_filename = self.handle_file_upload()

        if forward_data is not None and reverse_data is not None:
            try:
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                if params.get('use_smart_name', True) and fw_filename:
                    name_no_ext = fw_filename.rsplit('.', 1)[0]
                    parts = name_no_ext.split('_')
                    # ดึงส่วนที่ 3 (Index 2
                    sample_id = parts[2] if len(parts) >= 3 else name_no_ext
                    header_id = f"{sample_id}_{timestamp}"
                else:
                    header_id = f"Consensus_Sequence_{timestamp}"

                # Display raw sequence information
                self.display_raw_sequence_info(forward_data, reverse_data)

                # Perform quality trimming
                forward_trimmed, reverse_trimmed, reverse_trimmed_rc, reverse_qual_rc = self.perform_quality_trimming(
                    forward_data, reverse_data, params)

                # Generate consensus 
                alignment_result = self.generate_consensus(forward_trimmed, reverse_trimmed, reverse_trimmed_rc, reverse_qual_rc, params, header_id)

                if alignment_result:
                    # Run BLAST analysis if requested
                    blast_df = self.run_blast_analysis(alignment_result.consensus, params)

                    # Provide download options
                    self.provide_downloads(alignment_result, forward_trimmed, reverse_trimmed, reverse_trimmed_rc, params, blast_df, 
                                           len(forward_data.sequence), len(reverse_data.sequence), header_id)

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

