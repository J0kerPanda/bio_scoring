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
    In table representation seq1 is considered to be placed horizontally and
    seq2 is placed vertically
   */
  def apply(seq1: Sequence, seq2: Sequence): ProcessingResult = {
    val s1 = seq1.content
    val s2 = seq2.content

    val res = fillScoreMatrix(createScoreMatrix(s1.length, s2.length), s1, s2)
    res.foreach(e => println(e.mkString(" ")))

    val (adj1, adj2) = adjustSequences(res, s1, s2)
    println(adj1)
    println(adj2)
    ProcessingResult(0, "", "")
  }

  private def createScoreMatrix(strSize1: Int, strSize2: Int): ScoreMatrix = {
    val res: ScoreMatrix = Array.ofDim(strSize2 + 1, strSize1 + 1)
    res.indices.foreach(i => res(i)(0) = i * gapPenalty) //first column
    res(0).indices.drop(1).foreach(j => res(0)(j) = j * gapPenalty) //first row
    res
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  @tailrec
  private def fillScoreMatrix(scoreMatrix: ScoreMatrix, s1: String, s2: String, i1: Int = 1, i2: Int = 1): ScoreMatrix = {
    if ((i2 == s2.length) && (i1 == s1.length)) {
      scoreMatrix
    } else  {
      scoreMatrix.indices.drop(i2).foreach { i =>
        scoreMatrix(i)(i1) = computeScore(scoreMatrix, s1, s2, i1, i) //column i1
      }
      scoreMatrix(0).indices.drop(i1 + 1).foreach { j =>
        scoreMatrix(i2)(j) = computeScore(scoreMatrix, s1, s2, j, i2) //row i2
      }
      fillScoreMatrix(scoreMatrix, s1, s2, min(i1 + 1, s1.length), min(i2 + 1, s2.length))
    }
  }

  private def adjustSequences(scoreMatrix: ScoreMatrix, s1: String, s2: String): (String, String) = {
    import SequenceProcessor.gapSymbol
    val scoreMatrix = fillScoreMatrix(createScoreMatrix(s1.length, s2.length), s1, s2)
    val b1 = new StringBuilder()
    val b2 = new StringBuilder()

    def adjust(i1: Int, i2: Int): (Int, Int) = {
      val (u2: Int, u1: Int, f) = List(
        (i2 - 1, i1 - 1, () => { b1.append(s1(i1 - 1)); b2.append(s2(i2 - 1)) }),
        (i2 - 1, i1, () => { b1.append(gapSymbol); b2.append(s2(i2 - 1)) }),
        (i2, i1 - 1, () => { b1.append(s1(i1 - 1)); b2.append(gapSymbol) })
      )
        .map(ind => {
          val res = ind -> scoreMatrix(ind._1)(ind._2)
          println(res._1._1, res._1._2, res._2)
          res
        })
        .maxBy(_._2)
        ._1

      f()
      (u1, u2)
    }

    @tailrec
    def adjustRec(i1: Int, i2: Int): Unit = {
      if ((i1 != 1) || (i2 != 1)) {
        val (u1, u2) = adjust(i1, i2)
        adjustRec(math.max(u1, 1), math.max(u2, 1))
      }
    }

    adjustRec(s1.length, s2.length)
    adjust(1, 1)
    (b1.reverse.toString(), b2.reverse.toString())
  }

  /*
    Indices are greater than the current position in the corresponding string by 1
   */
  private def computeScore(scoreMatrix: ScoreMatrix, s1: String, s2: String, i1: Int, i2: Int): Int = {
    val diffScore = weightMatrix((s1(i1 - 1), s2(i2 - 1)))
    max(
      scoreMatrix(i2)(i1 - 1) + gapPenalty,
      scoreMatrix(i2 - 1)(i1) + gapPenalty,
      scoreMatrix(i2 - 1)(i1 - 1) + diffScore
    )
  }

  private def max(i1: Int, i: Int*): Int = math.max(i1, i.max)

  private def min(i1: Int, i: Int*): Int = math.min(i1, i.min)
}
