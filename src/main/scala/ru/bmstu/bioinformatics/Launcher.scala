package ru.bmstu.bioinformatics

import java.io.File

import ru.bmstu.bioinformatics.reader.fasta.FileReader
import ru.bmstu.bioinformatics.scoring.{SequenceProcessor, WeightMatrix}
import ru.bmstu.bioinformatics.sequence.{Nucleotide, Protein, Sequence}

object Launcher extends App {

  override def main(args: Array[String]): Unit = {

    val fastaFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/uniprot.fasta")
    val matrixFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/protein.mtx")

    val m = WeightMatrix(Nucleotide)

    val processor = new SequenceProcessor(-5, m)

//    println(processor(
//      Sequence("", "ATGAC"),
//      Sequence("", "AGGA")
//    ))

    println(processor(
      Sequence("", "ATTAAAT"),
      Sequence("", "AAAT")
    ))
  }
}
