package ru.bmstu.bioinformatics.scoring

import ru.bmstu.bioinformatics.scoring.WeightMatrix.KeyMatrix
import ru.bmstu.bioinformatics.sequence.Sequence
import SequenceProcessor._

import scala.annotation.tailrec

object SequenceProcessor {
  private val gapSymbol = '-'

  /** Class, representing leftmost column and upmost row of the matrix
    */
  trait MatrixCorner[A] {
    def vertSize: Int
    def horSize: Int

    def apply(i: Int)(j: Int): A
  }

  /** Corner element is present only in vertical vector
    */
  private case class MatrixCornerVertBiased[A] private(protected val vert: Vector[A],
                                                       protected val hor: Vector[A]) extends MatrixCorner[A] {

    override val vertSize: Int = vert.size
    override val horSize: Int = hor.size + 1

    override def apply(i: Int)(j: Int): A = if (j == 0) vert(i) else hor(j)
  }

  private object MatrixCornerVertBiased {

    def apply[A](fillRow: A, fillCol: A, rows: Int, cols: Int): MatrixCornerVertBiased[A] =
      MatrixCornerVertBiased(
        vert = Vector.fill(rows)(fillRow),
        hor = Vector.fill(cols - 1)(fillCol)
      )

    def apply[A](fillRow: Int => A, fillCol: Int => A, rows: Int, cols: Int): MatrixCornerVertBiased[A] =
      MatrixCornerVertBiased(
        vert = (0 until rows).map(fillRow).toVector,
        hor = (1 until cols).map(fillCol).toVector
      )
  }

  /** Corner element is present only in horizontal vector
    */
  private case class MatrixCornerHorBiased[A] private(protected val vert: Vector[A],
                                                      protected val hor: Vector[A]) extends MatrixCorner[A] {

    override val vertSize: Int = vert.size
    override val horSize: Int = hor.size + 1

    override def apply(i: Int)(j: Int): A = if (i == 0) hor(j) else vert(i)
  }

  private object MatrixCornerHorBiased {

    def apply[A](fillRow: A, fillCol: A, rows: Int, cols: Int): MatrixCornerVertBiased[A] =
      MatrixCornerVertBiased(
        vert = Vector.fill(rows - 1)(fillRow),
        hor = Vector.fill(cols)(fillCol)
      )

    def apply[A](fillRow: Int => A, fillCol: Int => A, rows: Int, cols: Int): MatrixCornerVertBiased[A] =
      MatrixCornerVertBiased(
        vert = (1 until rows).map(fillRow).toVector,
        hor = (0 until cols).map(fillCol).toVector
      )
  }

  private type PathPoint = (Int, Int)
  private type ScorePathCorner = MatrixCornerVertBiased[(Int, List[PathPoint])]

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
                        continuousGapPenalty: Either[None.type, Int] = Left(None),
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
      matching = MatrixCornerVertBiased(Integer.MAX_VALUE, Integer.MAX_VALUE, rows, cols),
      gapS1 = continuousGapPenalty
        .map(cgp => MatrixCornerVertBiased(i => gapPenalty + cgp * (i - 1), _ => Integer.MAX_VALUE, rows, cols))
        .getOrElse(MatrixCornerVertBiased(gapPenalty, Integer.MAX_VALUE, rows, cols)),
      gapS2 = continuousGapPenalty
        .map(cgp => MatrixCornerHorBiased(_ => Integer.MAX_VALUE, j => gapPenalty + cgp * (j - 1),  rows, cols))
        .getOrElse(MatrixCornerHorBiased(Integer.MAX_VALUE, gapPenalty, rows, cols)),
    )

    val acc = MatrixCornerVertBiased((0, List.empty[PathPoint]), (0, List.empty[PathPoint]), rows, cols)

    val (score, path) = computeScorePath(initialTriplet, acc, s1, s2)
    val (adj1, adj2) = adjustSequences(s1, s2, path)
    ProcessingResult(score, adj1, adj2)
  }

  //Path is reversed
  @tailrec
  private def computeScorePath(currentTriplet: ProcessorTriplet,
                               acc: ScorePathCorner,
                               seq1: String,
                               seq2: String): (Int, List[PathPoint]) = {
    val (newTriplet, newAcc) = computeStep(currentTriplet, acc, seq1, seq2)

    if ((acc.vert.size == 1) && (acc.hor.size == 1)) {
      newAcc.vert(0)
    } else {
      computeScorePath(newTriplet, newAcc, seq1, seq2)
    }
  }

  private def computeStep(currentTriplet: ProcessorTriplet,
                          acc: ScorePathCorner,
                          seq1: String,
                          seq2: String): (ProcessorTriplet, ScorePathCorner) = {

    val builderMatching = Array.
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
