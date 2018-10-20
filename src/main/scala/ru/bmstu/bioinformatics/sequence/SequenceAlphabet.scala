package ru.bmstu.bioinformatics.sequence

object SequenceAlphabet {

  private val nucleotideAlphabet = Set('A', 'C', 'G', 'T')
  private val proteinAlphabet = Set('A' to 'Z':_*)
//  (
//    'A', 'B', 'C', 'D', 'E', 'F',
//    'G', 'H', 'I', 'K', 'L', 'M',
//    'N', 'P', 'Q', 'R', 'S', 'T',
//    'U', 'V', 'W', 'Y', 'X', 'Z'
//  )

  def apply(sequenceType: SequenceType): Set[Char] = sequenceType match {
    case Nucleotide => nucleotideAlphabet
    case Protein => proteinAlphabet
  }
}
