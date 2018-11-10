package ru.bmstu.bioinformatics.scoring

import ru.bmstu.bioinformatics.scoring.WeightMatrix.KeyMatrix
import ru.bmstu.bioinformatics.sequence.Sequence
import SequenceProcessor._

import scala.annotation.tailrec

object SequenceProcessor {
  private val gapSymbol = '-'

  /** Class, representing leftmost column and upmost row of the matrix
    * Corner element is present only in vertical vector
    */
  private case class MatrixCorner[A](vert: Vector[A],
                                     hor: Vector[A]) {

    val vertSize: Int = vert.size
    val horSize: Int = hor.size + 1

    def apply(i: Int)(j: Int): A = if (j == 0) vert(i) else hor(j-1)

    def print(): Unit = {
      println((vert(0) +: hor).mkString(","))
      vert.drop(1).foreach(println)
    }
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
    * @param gapS1 - corner of matrix responsible for first sequence gaps
    * @param gapS2 - corner of matrix responsible for second sequence gaps
    */
  private case class ProcessorPair(gapS1: MatrixCorner[Int],
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
                        continuousGapPenalty: Int = 0) {

  val limitValue: Int = Int.MinValue - 2 * gapPenalty + 1

  /*
    In corresponding table representation seq1 is considered to be placed vertically and
    seq2 is placed horizontally
   */
  def apply(seq1: Sequence, seq2: Sequence): ProcessingResult = {
    val s1 = seq1.content
    val s2 = seq2.content
    val rows = s1.length + 1
    val cols = s2.length + 1

    val pair = ProcessorPair(
      gapS1 =
        MatrixCorner(i => gapPenalty + continuousGapPenalty * (i - 1), _ => limitValue, rows, cols),
      gapS2 =
        MatrixCorner.applyHorizontalBias(_ => limitValue, j => gapPenalty + continuousGapPenalty * (j - 1),  rows, cols)
    )

    val acc: ScorePathCorner = MatrixCorner(
      pair.gapS1.vert.updated(0, 0).zipWithIndex.map { case (s, i) => (s, List.fill(i)((1, 0))) },
      pair.gapS2.hor.zipWithIndex.map { case (s, j) => (s, List.fill(j + 1)((0, 1))) }
    )

    val (score, path) = computeScorePath(pair, acc, s1, s2)
    val (adj1, adj2) = adjustSequences(s1, s2, path)
    ProcessingResult(score, adj1, adj2)
  }

  //Path is reversed
  @tailrec
  private def computeScorePath(currentPair: ProcessorPair,
                               acc: ScorePathCorner,
                               seq1: String,
                               seq2: String): (Int, List[PathPoint]) = {
    val (newPair, newAcc) = computeStep(currentPair, acc, seq1, seq2)

    if (newAcc.vertSize == 2) {
      val (_, _, sp) = computeHorizontal(newPair, newAcc, seq1, seq2)
      sp.last
    } else if (newAcc.horSize == 2) {
      val (_, _, sp) = computeVertical(newPair, newAcc, seq1, seq2)
      sp.last
    } else {
      computeScorePath(newPair, newAcc, seq1, seq2)
    }
  }

  private def createInitialArrays(size: Int,
                                  prevPair: ProcessorPair,
                                  acc: ScorePathCorner,
                                  charS1: Char,
                                  charS2: Char): (Array[Int], Array[Int], Array[Int], Array[(Int, List[PathPoint])]) = {

    val ProcessorPair(prevS1G, prevS2G) = prevPair
    val base = Array.fill(size)(limitValue)
    //Matching, sequence 1 gaps, sequence 2 gaps, score path
    val (m, s1g, s2g, sp) = (
      base.clone(),
      base.clone(),
      base.clone(),
      Array.ofDim[(Int, List[PathPoint])](size)
    )

    //Diagonal element score
    //matrix index, score, path
    val (mtx, score, path) = List(
      (m, acc(0)(0)._1 + weightMatrix((charS1, charS2)), (1, 1) :: acc(0)(0)._2),
      (s1g, math.max(acc(0)(1)._1 + gapPenalty, prevS1G(0)(1) + continuousGapPenalty), (1, 0) :: acc(1)(0)._2),
      (s2g, math.max(acc(1)(0)._1 + gapPenalty, prevS2G(1)(0) + continuousGapPenalty), (0, 1) :: acc(0)(1)._2)
    ).maxBy(_._2)

    mtx.update(0, score)
    sp.update(0, score -> path)

    (m, s1g, s2g, sp)
  }

  //Computes new column (all elements)
  private def computeVertical(prevPair: ProcessorPair,
                              acc: ScorePathCorner,
                              seq1: String,
                              seq2: String): (Vector[Int], Vector[Int], Vector[(Int, List[PathPoint])]) = {

    val vertOffset = math.min(seq1.length - 1, seq1.length + 1 - acc.vertSize)
    val horOffset = math.min(seq2.length - 1, seq2.length - acc.horSize + 1)

    val newVertSize = math.max(1, acc.vertSize - 1)

    //Matching vertical, sequence 1 vertical gaps, sequence 2 vertical gaps, score path vertical
    val (mv, s1gv, s2gv, spv) = createInitialArrays(newVertSize, prevPair, acc, seq1(vertOffset), seq2(horOffset))

    //process each row
    for {
      i <- 1 until newVertSize
    } {
      val (mtx, score, path) = List(
        (mv, acc(i)(0)._1 + weightMatrix((seq1(vertOffset + i), seq2(horOffset))), (1, 1) :: acc(i)(0)._2),
        (s1gv, math.max(spv(i - 1)._1 + gapPenalty, s1gv(i - 1) + continuousGapPenalty), (1, 0) :: spv(i - 1)._2),
        (s2gv, math.max(acc(i + 1)(0)._1 + gapPenalty, prevPair.gapS2(i + 1)(0) + continuousGapPenalty), (0, 1) :: acc(i + 1)(0)._2)
      ).maxBy(_._2)

      mtx.update(i, score)
      spv.update(i, score -> path)
    }

    (s1gv.toVector, s2gv.toVector, spv.toVector)
  }

  //Computes new row (all elements)
  private def computeHorizontal(prevPair: ProcessorPair,
                                acc: ScorePathCorner,
                                seq1: String,
                                seq2: String): (Vector[Int], Vector[Int], Vector[(Int, List[PathPoint])]) = {

    val vertOffset = math.min(seq1.length - 1, seq1.length + 1 - acc.vertSize)
    val horOffset = math.min(seq2.length - 1, seq2.length - acc.horSize + 1)

    val newHorSize = math.max(1, acc.horSize - 1)

    //Matching horizontal, sequence 1 horizontal gaps, sequence 2 horizontal gaps, score path horizontal
    val (mh, s1gh, s2gh, sph) = createInitialArrays(newHorSize, prevPair, acc, seq1(vertOffset), seq2(horOffset))
    //process each column

    for {
      j <- 1 until newHorSize
    } {
      val (mtx, score, path) = List(
        (mh, acc(0)(j)._1 + weightMatrix((seq1(vertOffset), seq2(horOffset + j))), (1, 1) :: acc(0)(j)._2),
        (s1gh, math.max(acc(0)(j + 1)._1 + gapPenalty, prevPair.gapS1(0)(j + 1) + continuousGapPenalty), (1, 0) :: acc(0)(j + 1)._2),
        (s2gh, math.max(sph(j - 1)._1 + gapPenalty, s2gh(j - 1) + continuousGapPenalty), (0, 1) :: sph(j - 1)._2)
      ).maxBy(_._2)

      mtx.update(j, score)
      sph.update(j, score -> path)
    }

    (s1gh.toVector, s2gh.toVector, sph.toVector)
  }

  private def computeStep(prevPair: ProcessorPair,
                          acc: ScorePathCorner,
                          seq1: String,
                          seq2: String): (ProcessorPair, ScorePathCorner) = {

    val (s1gv, s2gv, spv) = computeVertical(prevPair, acc, seq1, seq2)
    val (s1gh, s2gh, sph) = computeHorizontal(prevPair, acc, seq1, seq2)

    val newTriplet = ProcessorPair(
      gapS1 = MatrixCorner(s1gv, s1gh.drop(1)),
      gapS2 = MatrixCorner(s2gv, s2gh.drop(1))
    )

    val newAcc = MatrixCorner(spv, sph.drop(1))

    (newTriplet, newAcc)
  }

  // Path is reversed
  private def adjustSequences(s1: String, s2: String, path: List[PathPoint]): (String, String) = {
    import SequenceProcessor.gapSymbol
    val b1 = new StringBuilder()
    val b2 = new StringBuilder()

    var s1i = 0
    var s2i = 0

    println(s1.length, s2.length, path.length)

    path.reverse.foreach {
      case (1, 1) =>
        b1.append(s1(s1i))
        b2.append(s2(s2i))
        s1i += 1
        s2i += 1

      case (1, 0) =>
        b1.append(s1(s1i))
        b2.append(gapSymbol)
        s1i += 1

      case (0, 1) =>
        b1.append(gapSymbol)
        b2.append(s2(s2i))
        s2i += 1
    }

    (b1.toString(), b2.toString())
  }
}
