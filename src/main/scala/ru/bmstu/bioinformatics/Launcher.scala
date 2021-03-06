package ru.bmstu.bioinformatics

import ru.bmstu.bioinformatics.cli.CliParser
import ru.bmstu.bioinformatics.reader.fasta.FileReader
import ru.bmstu.bioinformatics.scoring.{SequenceProcessor, WeightMatrix}

import scala.util.{Failure, Success, Try}

object Launcher extends App {

  override def main(args: Array[String]): Unit = {

    CliParser(args) match {
      case Some(config) =>
        val res = for {
          s1 <- Try(FileReader.readFirst(config.seq1, config.sequenceType))
          s2 <- Try(FileReader.readFirst(config.seq2, config.sequenceType))
        } yield {
          val weightMatrix = WeightMatrix(config.sequenceType)
          new SequenceProcessor(config.gapPenalty, weightMatrix, config.sequenceGapPenalty.getOrElse(config.gapPenalty))(s1, s2)
        }

        res match {
          case Failure(e) =>
            println(e.getMessage)
            e.printStackTrace()

          case Success(result) =>
            result.print(config.groupSize)
        }

      case None =>
    }
  }
}
