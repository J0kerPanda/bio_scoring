package ru.bmstu.bioinformatics.scoring

import ru.bmstu.bioinformatics.scoring.SequenceProcessor.ProcessingResult
import ru.bmstu.bioinformatics.scoring.WeightMatrix.KeyMatrix
import ru.bmstu.bioinformatics.sequence.Sequence

import scala.annotation.tailrec

object SequenceProcessor {
  private val gapSymbol = '-'

  case class ProcessingResult(score: Int, adjustedSeq1: String, adjustedSeq2: String) {

    def print(groupSize: Int): Unit = {
      println(s"Score: $score")
      adjustedSeq1.grouped(groupSize)
        .zip(adjustedSeq2.grouped(groupSize))
        .zipWithIndex
        .foreach { case ((s1, s2), i) =>
          println(s"${i * groupSize + 1} - ${i * groupSize + s1.length}")
          println(s"s1: $s1")
          println(s"s2: $s2")
        }
    }
  }
}

class SequenceProcessor(gapPenalty: Int, weightMatrix: KeyMatrix, continuousGapPenalty: Int) {

  private type PathPoint = (Int, Int)
  private type ScorePathEntry = (Int, List[PathPoint])
  private type ScorePathMatrix = Array[Array[ScorePathEntry]]

  /*
    In corresponding table representation seq1 is considered to be placed vertically and
    seq2 is placed horizontally
   */
  def apply(seq1: Sequence, seq2: Sequence): ProcessingResult = {
    val s1 = seq1.content
    val s2 = seq2.content

    val scoreMatrix = fillScorePathMatrix(createScoreMatrix(s1.length, s2.length), s1, s2)
    val res = scoreMatrix.last.last
    val (adj1, adj2) = adjustSequences(s1, s2, res._2)
    ProcessingResult(res._1, adj1, adj2)
  }

  private def createScoreMatrix(strSize1: Int, strSize2: Int): ScorePathMatrix = {
    val res: ScorePathMatrix = Array.ofDim(strSize1 + 1, strSize2 + 1)
    res.indices.foreach(i => res(i)(0) = (i * gapPenalty, Nil)) //first column
    res(0).indices.drop(1).foreach(j => res(0)(j) = (j * gapPenalty, Nil)) //first row
    res
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  @tailrec
  private def fillScorePathMatrix(matrix: ScorePathMatrix, s1: String, s2: String, i1: Int = 1, i2: Int = 1): ScorePathMatrix = {
    matrix.indices.drop(i1).foreach { i =>
      matrix(i)(i2) = computeScorePath(matrix, s1, s2, i, i2) //column i2
    }
    matrix(0).indices.drop(i2 + 1).foreach { j =>
      matrix(i1)(j) = computeScorePath(matrix, s1, s2, i1, j) //row i1
    }

    if ((i1 != s1.length) || (i2 != s2.length)) {
      fillScorePathMatrix(matrix, s1, s2, math.min(i1 + 1, s1.length), math.min(i2 + 1, s2.length))
    } else {
      matrix
    }
  }

  // Path is reversed
  private def adjustSequences(s1: String, s2: String, path: List[PathPoint]): (String, String) = {
    import SequenceProcessor.gapSymbol
    val b1 = new StringBuilder()
    val b2 = new StringBuilder()

    val shiftedPath = (s1.length, s2.length) :: path.dropRight(1)

    shiftedPath.zip(path).foreach {
      case ((c1, c2), (o1, o2)) => (c1 - o1, c2 - o2) match {
        case (1, 1) =>
          b1.append(s1(o1))
          b2.append(s2(o2))

        case (1, 0) =>
          b1.append(s1(o1))
          b2.append(gapSymbol)

        case (0, 1) =>
          b1.append(gapSymbol)
          b2.append(s2(o2))
      }
    }


    (b1.reverse.toString(), b2.reverse.toString())
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  private def computeScorePath(matrix: ScorePathMatrix, s1: String, s2: String, i1: Int, i2: Int): ScorePathEntry = {
    val diffScore = weightMatrix((s1(i1 - 1), s2(i2 - 1)))
    val (pathPoint, score) = List(
      (i1, i2 - 1) -> (matrix(i1)(i2 - 1)._1 + gapPenalty),
      (i1 - 1, i2) -> (matrix(i1 - 1)(i2)._1 + gapPenalty),
      (i1 - 1, i2 - 1) -> (matrix(i1 - 1)(i2 - 1)._1 + diffScore)
    )
      .maxBy(_._2)

    // Path is added backwards
    (score, pathPoint :: matrix(pathPoint._1)(pathPoint._2)._2)
  }
}
