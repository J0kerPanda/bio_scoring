package ru.bmstu.bioinformatics.sequence

object SequenceAlphabet {

  private val nucleotideAlphabet = Set('A', 'C', 'G', 'T')
  private val proteinAlphabet = Set('A' to 'Z':_*)

  def apply(sequenceType: SequenceType): Set[Char] = sequenceType match {
    case Nucleotide => nucleotideAlphabet
    case Protein => proteinAlphabet
  }
}
