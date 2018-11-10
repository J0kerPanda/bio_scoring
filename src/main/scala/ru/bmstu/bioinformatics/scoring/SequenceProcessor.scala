package ru.bmstu.bioinformatics.scoring

import ru.bmstu.bioinformatics.scoring.WeightMatrix.KeyMatrix
import ru.bmstu.bioinformatics.sequence.Sequence
import SequenceProcessor._

import scala.annotation.tailrec

object SequenceProcessor {
  private val gapSymbol = '-'

  /** Class, representing leftmost column and upmost row of the matrix
    * Corner element is present only in down vector
    */
  private case class MatrixCorner[A](row: Vector[A], col: Vector[A])

  private type PathPoint = (Int, Int)

  private object MatrixCorner {

    def apply[A](fill: A, rows: Int, cols: Int): MatrixCorner[A] =
      MatrixCorner(Vector.fill(rows)(fill), Vector.fill(cols - 1)(fill))
  }

  /** Representation of intermediate step of score computation process
    * @param matching - corner of matrix responsible for matching
    * @param gapS1 - corner of matrix responsible for first sequence gaps
    * @param gapS2 - corner of matrix responsible for second sequence gaps
    */
  private case class ProcessorTriplet(matching: MatrixCorner[Int],
                                      gapS1: MatrixCorner[Int],
                                      gapS2: MatrixCorner[Int])

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

class SequenceProcessor(gapPenalty: Int,
                        weightMatrix: KeyMatrix,
                        continuousGapPenalty: Int = 0,
                        considerStartGaps: Boolean = false,
                        considerEndGaps: Boolean = false) {

  /*
    In corresponding table representation seq1 is considered to be placed vertically and
    seq2 is placed horizontally
   */
  def apply(seq1: Sequence, seq2: Sequence): ProcessingResult = {
    val s1 = seq1.content
    val s2 = seq2.content
    val rows = s1.length
    val cols = s2.length

    val initialTriplet = ProcessorTriplet(
      matching = MatrixCorner(0, rows, cols),
      gapS1 = MatrixCorner(0, rows, cols),
      gapS2 = MatrixCorner(0, rows, cols)
    )

    val acc = MatrixCorner((0, List.empty[PathPoint]), rows, cols)

    val (score, path) = computeScorePath(initialTriplet, acc, s1, s2)
    val (adj1, adj2) = adjustSequences(s1, s2, path)
    ProcessingResult(score, adj1, adj2)
  }

  //Path is reversed
  @tailrec
  private def computeScorePath(currentTriplet: ProcessorTriplet,
                               acc: MatrixCorner[(Int, List[PathPoint])],
                               seq1: String,
                               seq2: String): (Int, List[PathPoint]) = {
    //todo compute new scorepath corner as well as triplet
    val newTriplet: ProcessorTriplet = ???
    val newAcc: MatrixCorner[(Int, List[PathPoint])] = ???

    if ((acc.col.size == 1) && (acc.row.size == 1)) {
      newAcc.col(0)
    } else {
      computeScorePath(newTriplet, newAcc, seq1, seq2)
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
}
