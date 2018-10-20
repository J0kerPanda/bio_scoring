package ru.bmstu.bioinformatics.reader.fasta

import java.io.{File, FileNotFoundException}

import scala.io.Source

class FastaReader {

  def apply(fileName: String): Iterator[FastaSequence] = {
    val file = new File(fileName)
    if (!file.exists() || !file.isFile) {
    } else {
      Source.fromFile(file)
        .getLines()
        .grouped(2)
        .
    }
  }

  def validateFile(file: File): Boolean = {
    if (file.exists() || file.isFile) {

    } else {
      throw new FileNotFoundException("sequence file not found")
    }
  }
}
