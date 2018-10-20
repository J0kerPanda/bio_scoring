package ru.bmstu.bioinformatics.scoring

import ru.bmstu.bioinformatics.scoring.SequenceProcessor.ProcessingResult
import ru.bmstu.bioinformatics.scoring.WeightMatrix.KeyMatrix
import ru.bmstu.bioinformatics.sequence.Sequence

object SequenceProcessor {

  case class ProcessingResult(score: Int, adjustedSeq1: String, adjustedSeq2: String)
}

class SequenceProcessor(gapPenalty: Int, weightMatrix: KeyMatrix) {

  private type ScoreMatrix = Array[Array[Int]]

  /*
    In table representation seq1 is considered to be placed horizontally and
    seq2 is placed vertically
   */
  def apply(seq1: Sequence, seq2: Sequence): ProcessingResult = {
    val s1 = seq1.content
    val s2 = seq2.content

    val res = process(s1, s2, 0, 0, createScoreMatrix(s1.length, s2.length))
    res.foreach(e => println(e.mkString(" ")))
    ProcessingResult(0, "", "")
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  private def process(s1: String, s2: String, i1: Int, i2: Int, scoreMatrix: ScoreMatrix): ScoreMatrix = {
    if ((i2 <= s2.length) && (i1 <= s1.length)) {
      scoreMatrix.indices.drop(i2).foreach { i =>
        scoreMatrix(i)(i1) = computeScore(scoreMatrix, s1, s2, i1, i) //column i1
      }
      scoreMatrix(0).indices.drop(i1 + 1).foreach { j =>
        scoreMatrix(i2)(j) = computeScore(scoreMatrix, s1, s2, j, i2) //row i2
      }
      process(s1, s2, max(i1, s1.length), max(i2, s2.length), scoreMatrix)
    } else {
      scoreMatrix
    }
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  private def computeScore(scoreMatrix: ScoreMatrix, s1: String, s2: String, i1: Int, i2: Int): Int = {
    val diffScore = weightMatrix((s1(i1 - 1), s2(i2 - 2)))
    max(
      scoreMatrix(i2)(i1 - 1) + gapPenalty,
      scoreMatrix(i2 - 1)(i1) + gapPenalty,
      scoreMatrix(i2 - 1)(i1 - 1) + diffScore
    )
  }

  private def createScoreMatrix(strSize1: Int, strSize2: Int): ScoreMatrix = {
    val res: ScoreMatrix = Array.ofDim(strSize2 + 1, strSize1 + 1)
    res.indices.foreach(i => res(i)(0) = i * gapPenalty) //first column
    res(0).indices.drop(2).foreach(j => res(0)(j) = j * gapPenalty) //first row
    res
  }

  private def max(i1: Int, i: Int*): Int = math.max(i1, i.max)
}
