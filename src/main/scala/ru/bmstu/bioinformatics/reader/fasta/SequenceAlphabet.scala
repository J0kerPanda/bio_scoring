package ru.bmstu.bioinformatics.reader.fasta

object SequenceAlphabet {

  private val nucleotideAlphabet = Set('A', 'C', 'G', 'T')
  private val proteinAlphabet = Set('A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V')

  def apply(sequenceType: SequenceType): Set[Char] = sequenceType match {
    case Nucleotide => nucleotideAlphabet
    case Protein => proteinAlphabet
  }
}
