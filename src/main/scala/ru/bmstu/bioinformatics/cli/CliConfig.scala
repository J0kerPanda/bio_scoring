package ru.bmstu.bioinformatics.cli

import java.io.File

import ru.bmstu.bioinformatics.sequence.{Protein, SequenceType}

case class CliConfig(seq1: File = new File(""),
                     seq2: File = new File(""),
                     sequenceType: SequenceType = Protein,
                     output: Option[File] = None,
                     groupSize: Int = 100,
                     gapPenalty: Int = -10,
                     sequenceGapPenalty: Option[Int] = None)
