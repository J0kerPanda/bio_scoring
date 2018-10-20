package ru.bmstu.bioinformatics.cli

import java.io.File

import ru.bmstu.bioinformatics.sequence.SequenceType

case class CliConfig(seq1: File,
                     seq2: File,
                     sequenceType: SequenceType,
                     output: Option[File] = None,
                     groupSize: Int = 100,
                     gapPenalty: Int = 10)
