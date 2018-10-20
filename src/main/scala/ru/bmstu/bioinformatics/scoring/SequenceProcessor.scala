package ru.bmstu.bioinformatics.scoring

import ru.bmstu.bioinformatics.scoring.SequenceProcessor.ProcessingResult
import ru.bmstu.bioinformatics.scoring.WeightMatrix.KeyMatrix
import ru.bmstu.bioinformatics.sequence.Sequence

import scala.annotation.tailrec

object SequenceProcessor {
  private val gapSymbol = '-'

  case class ProcessingResult(score: Int, adjustedSeq1: String, adjustedSeq2: String)
}

class SequenceProcessor(gapPenalty: Int, weightMatrix: KeyMatrix) {

  private type ScoreMatrix = Array[Array[Int]]

  /*
    In corresponding table representation seq1 is considered to be placed vertically and
    seq2 is placed horizontally
   */
  def apply(seq1: Sequence, seq2: Sequence): ProcessingResult = {
    val s1 = seq1.content
    val s2 = seq2.content

    val scoreMatrix = fillScoreMatrix(createScoreMatrix(s1.length, s2.length), s1, s2)
    val (adj1, adj2) = adjustSequences(scoreMatrix, s1, s2)
    ProcessingResult(scoreMatrix.last.last, adj1, adj2)
  }

  private def createScoreMatrix(strSize1: Int, strSize2: Int): ScoreMatrix = {
    val res: ScoreMatrix = Array.ofDim(strSize1 + 1, strSize2 + 1)
    res.indices.foreach(i => res(i)(0) = i * gapPenalty) //first column
    res(0).indices.drop(1).foreach(j => res(0)(j) = j * gapPenalty) //first row
    res
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  @tailrec
  private def fillScoreMatrix(scoreMatrix: ScoreMatrix, s1: String, s2: String, i1: Int = 1, i2: Int = 1): ScoreMatrix = {
    if ((i1 == s1.length) && (i2 == s2.length)) {
      scoreMatrix
    } else  {
      scoreMatrix.indices.drop(i1).foreach { i =>
        scoreMatrix(i)(i2) = computeScore(scoreMatrix, s1, s2, i, i2) //column i2
      }
      scoreMatrix(0).indices.drop(i2 + 1).foreach { j =>
        scoreMatrix(i1)(j) = computeScore(scoreMatrix, s1, s2, i1, j) //row i1
      }
      fillScoreMatrix(scoreMatrix, s1, s2, min(i1 + 1, s1.length), min(i2 + 1, s2.length))
    }
  }

  private def adjustSequences(scoreMatrix: ScoreMatrix, s1: String, s2: String): (String, String) = {
    import SequenceProcessor.gapSymbol
    val scoreMatrix = fillScoreMatrix(createScoreMatrix(s1.length, s2.length), s1, s2)
    val b1 = new StringBuilder()
    val b2 = new StringBuilder()

    @tailrec
    def adjustRec(i1: Int = 0, i2: Int = 0): Unit = {
      val (u1: Int, u2: Int, f) = List(
        (i1 + 1, i2 + 1, () => { b1.append(s1(i1)); b2.append(s2(i2)) }),
        (i1 + 1, i2, () => { b1.append(s1(i1)); b2.append(gapSymbol) }),
        (i1, i2 + 1, () => { b1.append(gapSymbol); b2.append(s2(i2)) })
      )
        .map(ind => ind -> scoreMatrix(ind._1)(ind._2))
        .maxBy(_._2)
        ._1

      f()
      (u1, u2) match {
        //str 1 exhausted
        case (_, _) if u1 == s1.length  =>
          (s1.length until s2.length).foreach { i =>
            b1.append(gapSymbol)
            b2.append(s2(i))
          }

        //str 2 exhausted
        case (_, _) if u2 == s2.length =>
          (s2.length until s1.length).foreach { i =>
            b1.append(s1(i))
            b2.append(gapSymbol)
          }

        case (_, _) =>
          adjustRec(u1, u2)
      }
    }

    adjustRec()
    (b1.toString(), b2.toString())
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  private def computeScore(scoreMatrix: ScoreMatrix, s1: String, s2: String, i1: Int, i2: Int): Int = {
    val diffScore = weightMatrix((s1(i1 - 1), s2(i2 - 1)))
    max(
      scoreMatrix(i1)(i2 - 1) + gapPenalty,
      scoreMatrix(i1 - 1)(i2) + gapPenalty,
      scoreMatrix(i1 - 1)(i2 - 1) + diffScore
    )
  }

  private def max(i1: Int, i: Int*): Int = math.max(i1, i.max)

  private def min(i1: Int, i: Int*): Int = math.min(i1, i.min)
}
