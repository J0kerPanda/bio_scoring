package ru.bmstu.bioinformatics

import java.io.File

import ru.bmstu.bioinformatics.scoring.{SequenceProcessor, WeightMatrix}
import ru.bmstu.bioinformatics.sequence.{Nucleotide, Protein, Sequence}

object Launcher extends App {

  override def main(args: Array[String]): Unit = {

    val fastaFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/uniprot.fasta")
    val matrixFile = new File("/home/antony/temprojects/bio_scoring/src/main/resources/protein.mtx")

    val m = WeightMatrix(Nucleotide)
    val m2 = WeightMatrix(Protein)

    val processor = new SequenceProcessor(-5, m)
    val processor2 = new SequenceProcessor(-5, m2)

    println(processor(
      Sequence("", "ATTA"),
      Sequence("", "ATTA")
    ))

    println(processor(
      Sequence("", "ATTAAAT"),
      Sequence("", "AAAT")
    ))

    println(processor(
      Sequence("", "AGTAAA"),
      Sequence("", "AGTA")
    ))

    processor2(
      Sequence("", "MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH"),
      Sequence("", "MVHLTDAEKAAVNGLWGKVNPDDVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPKVKAHGKKVINAFNDGLKHLDNLKGTFAHLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGKEFTPCAQAAFQKVVAGVASALAHKYH"
    )).print(50)
  }
}
