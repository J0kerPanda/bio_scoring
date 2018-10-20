package ru.bmstu.bioinformatics.reader.fasta

sealed trait SequenceType

case object Nucleotide
case object Protein
