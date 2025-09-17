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
- Quality Cutoff: 20 (Minimum phred score threshold) (Sequence cannot be trimmed by more than 30% of the sequence length)\
- Window Size: 10 (Size of sliding window for quality assessment)
- Adjustment: 0 (How much bp further must the Best Window be to override the First Acceptable Window? If the adjustment is 0, select the First Acceptable Window.)
  
  <img width="465" height="567" alt="image" src="https://github.com/user-attachments/assets/8c1be7a2-62b3-4bea-bdd3-a3c47bd7ad26" />

### Results
- 🔁 Reverse complement + reverse Phred score handling
  <img width="1387" height="800" alt="image" src="https://github.com/user-attachments/assets/81bc6e28-daff-4a9f-9dee-1d635bcda4e9" />
- 🧩 Overlap alignment
  <img width="1452" height="674" alt="image" src="https://github.com/user-attachments/assets/7935d5fc-ef6d-4e0a-875a-328ae631ff57" />
- 📊 Output stats
  |   Metric  | Value |
  | --------- | ----- |
  |Matches |  277|
  |Mismatches|  12|
  |Gaps|  1101|
  |Alignment Score|  227.5|
  |Consensus_Length|  1390|
  |Quality_Cutoff|  20|
  |Window_Size|  10|
  |Adjustment|  0|

- 💾 Export to FASTA, CSV
  <img width="1186" height="179" alt="image" src="https://github.com/user-attachments/assets/228b01c5-1571-40b4-b679-dde72112508d" />
- 💾 (Optional) Blast Result
  <img width="1187" height="747" alt="image" src="https://github.com/user-attachments/assets/19abcc7b-470e-42cc-b1f1-dd57df1ad37c" />
- 🌐 Web interface (Streamlit) — runs locally or on Streamlit Cloud

---


