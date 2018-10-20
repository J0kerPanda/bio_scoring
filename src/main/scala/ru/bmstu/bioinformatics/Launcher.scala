package ru.bmstu.bioinformatics

import java.io.File

import ru.bmstu.bioinformatics.reader.fasta.FileReader
import ru.bmstu.bioinformatics.scoring.WeightMatrix
import ru.bmstu.bioinformatics.sequence.Protein

object Launcher extends App {

  override def main(args: Array[String]): Unit = {

    val fastaFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/uniprot.fasta")
    val matrixFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/matrix.txt")

    val m = WeightMatrix.fromFile(matrixFile)
    m.foreach(println)
    println(m.size)
  }
}
