# 🔬 SangerApp

A lightweight, Python-based Streamlit web application for generating high-quality consensus sequences from paired-end Sanger `.ab1` chromatogram files.

Designed for use in molecular biology, diagnostics, and teaching labs.

Link for running program: https://sangerapp.streamlit.app/

---

## ✨ Tutorial

### 🧬 Input paired `.ab1` Sanger reads (forward and reverse)
  <img width="1225" height="603" alt="image" src="https://github.com/user-attachments/assets/b0f24f27-62cc-409e-a171-4df269c2c16b" />

### ✂️ Automatic quality-based trimming
  
#### Analysis Parameters Default:
- Quality Cutoff: 20 (Minimum phred score threshold) (Sequence cannot be trimmed by more than 30% of the sequence length)
- Window Size: 10 (Size of sliding window for quality assessment)
- Adjustment: 0 (How much bp further must the Best Window be to override the First Acceptable Window? If the adjustment is 0, select the First Acceptable Window.)
  
  <img width="465" height="567" alt="image" src="https://github.com/user-attachments/assets/8c1be7a2-62b3-4bea-bdd3-a3c47bd7ad26" />

### Results
- 🔁 Reverse complement + reverse Phred score handling
  <img width="1512" height="490" alt="image" src="https://github.com/user-attachments/assets/984ca118-234c-4f47-8979-7fbf49825a1e" />
- 🧩 Overlap alignment
  <img width="1515" height="657" alt="image" src="https://github.com/user-attachments/assets/bf2b5900-1b3c-4d7f-8565-558abccd817c" />
- 📊 Output stats
  |   Metric  | Value |
  | --------- | ----- |
  |Matches |  423|
  |Mismatches|  0|
  |Gaps|  1005|
  |Alignment Score|  408|
  |Consensus_Length|  1428|
  |Quality_Cutoff|  20|
  |Window_Size|  10|
  |Start_Adjustment_5_end|  0|
  |Threshold_Type|  Window Mean|
  |Forward_Raw_Length|  1049|
  |Forward_Trimmed_Length|  937|
  |Forward_Kept_Region_bp|  19 - 955|
  |Reverse_Raw_Length|  1009|
  |Reverse_Trimmed_Length|  914|
  |Reverse_Kept_Region_bp|  24 - 937|

- 💾 Export to FASTA, CSV
  <img width="1186" height="179" alt="image" src="https://github.com/user-attachments/assets/228b01c5-1571-40b4-b679-dde72112508d" />
- 💾 (Optional) Blast Result
  <img width="1487" height="582" alt="image" src="https://github.com/user-attachments/assets/a7fd8034-7207-4a66-b30c-e377c1ecd1d2" />
- 🌐 Web interface (Streamlit) — runs locally or on Streamlit Cloud

---


