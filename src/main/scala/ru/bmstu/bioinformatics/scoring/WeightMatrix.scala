package ru.bmstu.bioinformatics.scoring

import java.io.{File, FileNotFoundException}

import scala.collection.mutable
import scala.io.Source

object WeightMatrix {

  type KeyMatrix = Map[(String, String), Int]

  private val proteinMatrix = Map(

  )

  /* Any lines starting with # are discarded
    The first line (excluding comments) should contain only space separated latin characters or
    an asterisk (any not-stated latin character)
    The body of the matrix (should be square) consists of lines starting with a latin character or
    an asterisk which is followed by a series of numbers corresponding to the size of the matrix
   */
  ////todo error handling
  def fromFile(file: File): KeyMatrix = {
    if (file.exists() && file.isFile) {
      val builder = mutable.HashMap()
      val names :: body = Source
        .fromFile(file)
        .getLines()
        .map(_.trim)
        .filter(_.startsWith("#"))
        .toList

      val nameArray = names.toCharArray.filter(c => (c != '*') && (c != ' '))
      val excludedSet = ('A' to 'Z').toSet.diff(nameArray.toSet)

      body.foreach { line =>
        val lineArray = line.split(" ")
        assert(lineArray(0).length == 1)
        val name = lineArray(0).head
        val weightsInt = lineArray.drop(1).map(_.toInt)
      }

      Map()

    } else {
      throw new FileNotFoundException(s"Matrix file not found [${file.getAbsolutePath}]")
    }
  }
}
