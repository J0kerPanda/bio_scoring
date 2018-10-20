package ru.bmstu.bioinformatics.sequence

sealed trait SequenceType

case object Nucleotide extends SequenceType
case object Protein extends SequenceType
