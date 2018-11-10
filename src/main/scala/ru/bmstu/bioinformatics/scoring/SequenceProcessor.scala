package ru.bmstu.bioinformatics.scoring

import ru.bmstu.bioinformatics.scoring.WeightMatrix.KeyMatrix
import ru.bmstu.bioinformatics.sequence.Sequence
import SequenceProcessor._

import scala.annotation.tailrec

object SequenceProcessor {
  private val gapSymbol = '-'
  private val limitValue = Integer.MIN_VALUE

  /** Class, representing leftmost column and upmost row of the matrix
    * Corner element is present only in vertical vector
    */
  private case class MatrixCorner[A](protected val vert: Vector[A],
                                     protected val hor: Vector[A]) {

    val vertSize: Int = vert.size
    val horSize: Int = hor.size + 1

    def apply(i: Int)(j: Int): A = if (j == 0) vert(i) else hor(j)
  }

  private object MatrixCorner {

    def apply[A](fill: A, rows: Int, cols: Int): MatrixCorner[A] =
      MatrixCorner(
        vert = Vector.fill(rows)(fill),
        hor = Vector.fill(cols - 1)(fill)
      )

    def apply[A](fillRow: Int => A, fillCol: Int => A, rows: Int, cols: Int): MatrixCorner[A] =
      MatrixCorner(
        vert = (0 until rows).map(fillRow).toVector,
        hor = (1 until cols).map(fillCol).toVector
      )

    def applyHorizontalBias[A](fillRow: Int => A, fillCol: Int => A, rows: Int, cols: Int): MatrixCorner[A] =
      MatrixCorner(
        vert = (0 until rows).map(fillRow).toVector.updated(0, fillCol(0)),
        hor = (1 until cols).map(fillCol).toVector
      )
  }

  private type PathPoint = (Int, Int)
  private type ScorePathCorner = MatrixCorner[(Int, List[PathPoint])]

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
      matching = MatrixCorner(limitValue, rows, cols),
      gapS1 =
        MatrixCorner(i => gapPenalty + continuousGapPenalty * (i - 1), _ => limitValue, rows, cols),
      gapS2 =
        MatrixCorner.applyHorizontalBias(_ => limitValue, j => gapPenalty + continuousGapPenalty * (j - 1),  rows, cols)
    )

    val acc = MatrixCorner((0, List.empty[PathPoint]), rows, cols)

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

    if ((acc.vertSize == 1) && (acc.horSize == 1)) {
      newAcc(0)(0)
    } else {
      computeScorePath(newTriplet, newAcc, seq1, seq2)
    }
  }

  private def computeStep(prevTriplet: ProcessorTriplet,
                          acc: ScorePathCorner,
                          seq1: String,
                          seq2: String): (ProcessorTriplet, ScorePathCorner) = {

    val s1Char = seq1(seq1.length - acc.vertSize)
    val s2Char = seq2(seq2.length - acc.horSize)

    //process each row
    val baseVertical = Array.fill(math.max(1, acc.vertSize - 1))(limitValue)
    //Matching vertical, sequence 1 vertical gaps, sequence 2 vertical gaps, score path vertical
    //s2gv nee
    val (mv, s1gv, s2gv, spv) = (baseVertical.clone(), baseVertical.clone(), baseVertical.clone(), baseVertical.clone())

    def compute

    for {
      i <- 1 to acc.vertSize
    } {
      val prevScore = acc(i)(0)
      val mScore = acc(i - 1)(0) ._1 + weightMatrix((s1Char, s2Char))
      val gapS1Score = acc(i - 1)(0) + if (prevTriplet.gapS1(i - 1)(0) != limitValue) continuousGapPenalty else gapPenalty
    }

    //process each column


    (???, ???)
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
