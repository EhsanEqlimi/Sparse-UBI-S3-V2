# Mixing Matrix Identification using k-SCA and S3

This repository contains a novel method for identifying the mixing matrix in underdetermined blind source separation (UBSS) based on k-sparse component analysis (k-SCA). The primary goal of this method is to recover the mixing matrix A when only the mixed signals X are known, with the additional constraint that the source signals S have k-sparse columns.

## Problem Statement

In UBSS, we have the following problem:

- We know the mixed signals X.
- However, both the mixing matrix A and the source signals S are unknown.
- There are more sources than sensors, making the problem underdetermined.
- The source signals S have k-sparse columns, meaning there are only k non-zero elements in each column.

## Aim

The main objective of this method is to find the mixing matrix A, given the known mixed signals X and the k-sparse constraint on the source signals S.

## Idea: Subspace Selective Search (S3)

This method is based on the Subspace Selective Search (S3) algorithm, which is designed to tackle the underdetermined blind source separation problem by exploiting the k-sparse nature of the source signals.

## References

If you use or reference this method, please make sure to cite the following papers:

1. [Multiple Sparse Component Analysis Based on Subspace Selective Search Algorithm]([link_to_paper_1](https://ieeexplore.ieee.org/abstract/document/7146277)) by E. Eqlimi and B. Makkiabadi, presented at the 2015 23rd Iranian Conference on Electrical Engineering.

2. [An Efficient K-SCA Based Underdetermined Channel Identification Algorithm for Online Applications]([link_to_paper_2](https://ieeexplore.ieee.org/document/7362867)) by E. Eqlimi and B. Makkiabadi, presented at the 2015 23rd European Signal Processing Conference (EUSIPCO).

## Contact Information

(C) Ehsan Eqlimi and Dr. Bahador Makkibadi, Jan 2015
Medical Physics & Biomedical Engineering Department
Tehran University of Medical Sciences (TUMS), Tehran, Iran

For questions or inquiries, please contact:
- Ehsan Eqlimi: [ehsan.eqlimi@outlook.com](mailto:ehsan.eqlimi@outlook.com)
- Dr. Bahador Makkibadi: [email@example.com](mailto:email@example.com)
