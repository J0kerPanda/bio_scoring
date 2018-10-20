package ru.bmstu.bioinformatics.reader.fasta

import java.io.{File, FileNotFoundException}

import ru.bmstu.bioinformatics.reader.fasta.FileReader.InvalidFileFormatException

import scala.io.Source

object FileReader {

  case class InvalidFileFormatException(file: File) extends Throwable {
    override def getMessage: String = s"Invalid sequence format in file [$file]"
  }
}

class FileReader {

  def apply(fileName: String, sequenceType: SequenceType): Iterator[Sequence] = {
    val file = new File(fileName)
    if (validateFile(file, sequenceType)) {
      groupedIterator(file)
        .map { case Seq(name, content) =>
          Sequence(name, content)
        }
    } else {
      throw InvalidFileFormatException(file)
    }
  }

  def validateFile(file: File, sequenceType: SequenceType): Boolean = {
    if (file.exists() || file.isFile) {
      val alphabet = SequenceAlphabet(sequenceType)
      !groupedIterator(file)
        .exists {
          case Seq(a, b) if a.startsWith(">") && !b.contains(!alphabet(_))=>
            false

          case _ =>
            true
        }
    } else {
      throw new FileNotFoundException(s"Sequence file [${file.getAbsolutePath}] not found")
    }
  }

  private def groupedIterator(file: File): Iterator[String]#GroupedIterator[String] = {
    Source.fromFile(file).getLines().grouped(2)
  }
}
