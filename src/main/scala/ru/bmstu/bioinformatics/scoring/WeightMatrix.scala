package ru.bmstu.bioinformatics.scoring

import java.io.{File, FileNotFoundException}

import scala.collection.mutable
import scala.io.Source

object WeightMatrix {

  type KeyMatrix = Map[(Char, Char), Int]

  /* Any lines starting with # are discarded
    The first line (excluding comments) should contain only space separated latin characters or
    an asterisk (any not-stated latin character)
    The body of the (square) matrix should consist of lines starting with a latin character or
    an asterisk which is followed by a series of numbers corresponding to the size of the matrix
   */
  //todo error handling
  def fromFile(file: File): KeyMatrix = {
    if (file.exists() && file.isFile) {
      val builder = mutable.HashMap[(Char, Char), Int]()
      val names :: body = Source
        .fromFile(file)
        .getLines()
        .map(_.trim)
        .filter(!_.startsWith("#"))
        .toList

      val nameArray = splitBySpaces(names).map { n =>
        assert(n.length == 1)
        n.head
      }
      assert(nameArray.length == body.length)
      val excludedSet = ('A' to 'Z').toSet.diff(nameArray.toSet)

      body.foreach { line =>
        val lineArray = splitBySpaces(line)
        assert(lineArray(0).length == 1)
        val rowName = lineArray(0).head

        lineArray
          .drop(1)
          .map(_.toInt)
          .zipWithIndex
          .foreach { case (weight, i) =>
            (rowName, nameArray(i)) match {
              case ('*', '*') =>
                for {
                  n1 <- excludedSet
                  n2 <- excludedSet
                } { builder.update((n1, n2), weight) }

              case ('*', colName) =>
                excludedSet.foreach(n => builder.update((n, colName), weight))

              case (_, '*') =>
                excludedSet.foreach(n => builder.update((rowName, n), weight))

              case (_, colName) =>
                builder.update((rowName, colName) , weight)
            }
          }
      }

      builder.toMap

    } else {
      throw new FileNotFoundException(s"Matrix file not found [${file.getAbsolutePath}]")
    }
  }

  private def splitBySpaces(str: String) = str.split("\\s+")
}
