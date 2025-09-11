# 🔬 Sanger Assembly App

A lightweight, Python-based Streamlit web application for generating high-quality consensus sequences from paired-end Sanger `.ab1` chromatogram files.

Designed for use in molecular biology, diagnostics, and teaching labs.

---

## ✨ Features

- 🧬 Input paired `.ab1` Sanger reads (forward and reverse)
- <img width="1225" height="603" alt="image" src="https://github.com/user-attachments/assets/b0f24f27-62cc-409e-a171-4df269c2c16b" />
- ✂️ Automatic quality-based trimming
- 🔁 Reverse complement + reverse Phred score handling
- 🧩 Overlap alignment using Clustal-style scoring
- 🧠 Quality-aware consensus: trust forward in first half, reverse in second
- 📊 Output stats: identity, coverage, mismatches, alignment score
- 💾 Export to FASTA, CSV, and optional BLAST-ready format
- 🌐 Web interface (Streamlit) — runs locally or on Streamlit Cloud

---

## 🚀 Getting Started

### 📦 Install requirements

```bash
pip install -r requirements.txt
