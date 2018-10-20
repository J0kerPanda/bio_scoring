package ru.bmstu.bioinformatics.cli

import java.io.File

import ru.bmstu.bioinformatics.sequence.{Nucleotide, Protein}
import scopt.OptionParser

object CliParser {

  def apply() = {
    new OptionParser[CliConfig]("bio_scoring") {

      opt[File]("seq1")
        .required()
        .text("First input sequence file path")
        .action((file, c) => c.copy(seq1 = file))

      opt[File]("seq2")
        .required()
        .text("Second input sequence file path")
        .action((file, c) => c.copy(seq2 = file))

      opt[Char]('s', "sequenceType")
        .required()
        .text("Type of sequences (n - nucleotides, p - proteins)")
        .validate {
          case 'n' | 'p' => success
          case _ => failure("Sequence type must be either n - nucleotides or p - proteins")
        }
        .action((s, c) => s match {
          case 'n' => c.copy(sequenceType = Nucleotide)
          case 'p' => c.copy(sequenceType = Protein)
        })

      opt[Int]('g', "gap")
        .optional()
        .text("Gap penalty")
        .action((p, c) => c.copy(gapPenalty = p))

      opt[Int]("group")
        .optional()
        .text("Length of groups to cut sequences into")
        .action((g, c) => c.copy(groupSize = g))
    }
  }
}
