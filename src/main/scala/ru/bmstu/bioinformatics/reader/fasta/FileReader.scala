package ru.bmstu.bioinformatics.reader.fasta

import java.io.{File, FileNotFoundException}

import scala.io.Source

object FileReader {

  case class InvalidFileFormatException(file: File) extends Throwable {
    override def getMessage: String = s"Invalid sequence format in file [$file]"
  }

  def readFirst(fileName: String, sequenceType: SequenceType): Option[Sequence] = {
    val iterator = FileReader(fileName, sequenceType)
    if (iterator.hasNext) {
      Some(iterator.next())
    } else {
      None
    }
  }

  def readAll(fileName: String, sequenceType: SequenceType): List[Sequence] = {
    FileReader(fileName, sequenceType).toList
  }

  def apply(fileName: String, sequenceType: SequenceType): Iterator[Sequence] = {
    val file = new File(fileName)
    if (validateFile(file, sequenceType)) {
      groupedIterator(file)
        .map { case (name, content) =>
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
          case (name, content) if name.startsWith(">") && content.nonEmpty && content.forall(alphabet(_)) =>
            false

          case _ =>
            true
        }
    } else {
      throw new FileNotFoundException(s"Sequence file [${file.getAbsolutePath}] not found")
    }
  }

  private def groupedIterator(file: File): Iterator[(String, String)] = {
    val baseIterator = Source.fromFile(file).getLines()

    new Iterator[(String, String)] {
      private var previousCaption: String = _

      override def hasNext: Boolean = baseIterator.hasNext

      override def next(): (String, String) = {
        val caption = if (previousCaption != null) previousCaption else baseIterator.next()
        val builder = new StringBuilder()
        var currLine = ""
        while (baseIterator.hasNext && !currLine.startsWith(">")) {
          currLine = baseIterator.next()
          builder.append(currLine)
        }
        previousCaption = currLine
        (caption, builder.toString())
      }
    }
  }
}
