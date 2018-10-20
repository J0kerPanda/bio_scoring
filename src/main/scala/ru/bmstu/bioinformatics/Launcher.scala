package ru.bmstu.bioinformatics

import ru.bmstu.bioinformatics.reader.fasta.FileReader
import ru.bmstu.bioinformatics.sequence.Protein

object Launcher extends App {

  override def main(args: Array[String]): Unit = {

    val fileName = "/home/antony/temprojects/bio_scoring/src/main/resources/uniprot.fasta"

    FileReader(fileName = fileName, sequenceType = Protein)
      .size
  }
}
