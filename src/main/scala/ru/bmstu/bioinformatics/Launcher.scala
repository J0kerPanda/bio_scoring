package ru.bmstu.bioinformatics

import ru.bmstu.bioinformatics.cli.CliParser
import ru.bmstu.bioinformatics.reader.fasta.FileReader
import ru.bmstu.bioinformatics.scoring.{SequenceProcessor, WeightMatrix}
import ru.bmstu.bioinformatics.sequence.{Protein, Sequence}

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
          new SequenceProcessor(config.gapPenalty, weightMatrix)(s1, s2)
        }

        res match {
          case Failure(e) =>
            println(e.getMessage)

          case Success(result) =>
            result.print(config.groupSize)
        }

      case None =>
    }

//    val fastaFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/uniprot.fasta")
//    val matrixFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/protein.mtx")
//
//    val m = WeightMatrix(Nucleotide)
//    val m2 = WeightMatrix(Protein)
//
//    val processor = new SequenceProcessor(-5, m)
//    val processor2 = new SequenceProcessor(-5, m2)
//
//    println(processor(
//      Sequence("", "ATTA"),
//      Sequence("", "ATTA")
//    ))
//
//    println(processor(
//      Sequence("", "ATTAAAT"),
//      Sequence("", "AAAT")
//    ))
//
//    println(processor(
//      Sequence("", "AGTAAA"),
//      Sequence("", "AGTA")
//    ))
//
//    processor2(
//      Sequence("", "MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH"),
//      Sequence("", "MVHLTDAEKAAVNGLWGKVNPDDVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPKVKAHGKKVINAFNDGLKHLDNLKGTFAHLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGKEFTPCAQAAFQKVVAGVASALAHKYH"
//    )).print(50)
  }
}
