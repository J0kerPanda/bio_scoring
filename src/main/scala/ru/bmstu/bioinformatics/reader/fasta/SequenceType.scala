package ru.bmstu.bioinformatics.reader.fasta

sealed trait SequenceType

case object Nucleotide extends SequenceType
case object Protein extends SequenceType
