package abeellab.parrot.rflp

import java.io.PrintWriter
import net.sf.samtools.SAMFileReader
import java.io.File
import net.sf.jannot.utils.SequenceTools
import org.arabidopsis.ahocorasick.AhoCorasick
import scala.collection.JavaConversions._
import org.arabidopsis.ahocorasick.SearchResult
import be.abeel.util.CountMap
import scala.io.Source
import atk.util.TimeInterval
import abeellab.parrot.Main

object RFLPTyper extends Main {

  case class Config(val outputFile: String = null, val input: File = null, val threshold: Int = 20, val flankSize: Int = 50, val genomeLen: Int = 4411708, val alignmentFraction: Double = 0.9, val verbose: Boolean = false)
  private val default = new Config

  def revcomp(read: Array[Byte]) = {
    val out = Array.ofDim[Byte](read.length)
    for (i <- 0 until read.length) {
      out(read.length - i - 1) = SequenceTools.complement(read(i).toChar).toByte
    }
    out
  }
  /* Load spacers */
  val primer5 = "ACTCACCGGGGCGGTT" /* from http://www.biomedcentral.com/1471-2164/13/249 */
  val primer3 = "TGAACCGCCCCGGTGA" /* Eyeballed from GenomeView */

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar parrot.jar rflp-is6110") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Read data in (unaligned) bam file.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("Name of output file.")
      opt[Int]('f', "flankSize") action { (x, c) => c.copy(flankSize = x) } text ("Minimum flank length to be included. Recommended 1/2 of read length default=" + default.flankSize)
      opt[Int]('t', "threshold") action { (x, c) => c.copy(threshold = x) } text ("Threshold for considering flank real. default=" + default.threshold)
      opt[Int]("genomeSize") action { (x, c) => c.copy(genomeLen = x) } text ("[Recommended not to change] Genome size for estimating depth. default=" + default.genomeLen)
      opt[Int]("alignmentFraction") action { (x, c) => c.copy(alignmentFraction = x) } text ("[Recommended not to change] Fraction of seq. library aligning to genome. default=" + default.alignmentFraction)
      opt[Unit]("verbose") action { (x, c) => c.copy(verbose = true) } text ("Loads of debug output")

    }

    parser.parse(args, Config()) map { config =>
      doStuff(config)
    }

  }
  def doStuff(config: Config) {

    /* Prep index */
    val ac = new AhoCorasick();

    ac.add(primer5.getBytes(), "primer5")
    ac.add(revcomp(primer5.getBytes()), "primer5_rc")
    ac.add(primer3.getBytes(), "primer3")
    ac.add(revcomp(primer3.getBytes()), "primer3_rc")

    ac.prepare();

    val pw = if (config.outputFile != null) new PrintWriter(config.outputFile) else new PrintWriter(System.out)

    val flanksOutput = if (config.outputFile != null) new File(config.outputFile + ".flanks") else { val f = File.createTempFile("parakeet", ".tmp"); f.deleteOnExit(); f }
    val pwFlanks = new PrintWriter(flanksOutput);
    pwFlanks.println(generatorInfo)

    pw.println(generatorInfo)

    /* Connect to bam file*/
    val inputSam = new SAMFileReader(config.input);
    val cm = new CountMap[String]
    //    val count = 0;

    /* Iterate over bamfile */
    val it = inputSam.iterator()
    var progress: Int = 0
    var totalCoverage: Long = 0

    val time = System.currentTimeMillis();
    val nano = System.nanoTime()
    while (it.hasNext()) {
      val sr = it.next()

      totalCoverage += sr.getReadLength()

      val read = sr.getReadBases()

      val result = ac.search(read);

      for (a <- result) {
        val result = a.asInstanceOf[SearchResult]
        val li = result.getLastIndex()
        val r = result.getOutputs().head.asInstanceOf[String]
        if (r.startsWith("primer5")) {
          if (r.endsWith("rc")) {
            // still to do
            pwFlanks.println("- 5" + new String(revcomp(sr.getReadBases())).substring(sr.getReadLength() + primer5.length() - li))

          } else {
            pwFlanks.println("+ 5" + sr.getReadString().substring(result.getLastIndex()));
          }
        }
        if (r.startsWith("primer3")) {
          if (!r.endsWith("rc")) {
            // still to do
            pwFlanks.println("+ 3" + new String(revcomp(sr.getReadBases())).substring(sr.getReadLength() + primer3.length() - li))

          } else {
            pwFlanks.println("- 3" + sr.getReadString().substring(result.getLastIndex()));
          }
        }
        
        
        assume(result.getOutputs().size == 1)
        for (s <- result.getOutputs()) {

          cm.count(s.asInstanceOf[String])

        }
        //       
      }

      // 2 snapshot coverage for all spacers through-out bam file
      progress += 1
      if (progress % 100000 == 0) {
        val expectedCoverage = (totalCoverage * config.alignmentFraction) / config.genomeLen
        println(progress + "\t" + cm + "\t" + expectedCoverage + "\t" + ((System.nanoTime() - nano) / progress / 1.0e6) + " ms/read\t" + new TimeInterval(System.currentTimeMillis() - time))
        if (config.verbose)
          pw.println(progress + "\t" + cm + "\t" + expectedCoverage + "\t" + ((System.nanoTime() - nano) / progress / 1.0e6) + " ms/read\t" + new TimeInterval(System.currentTimeMillis() - time))

      }
    }

    /* close bam file */
    it.close
    /* report results */
    if (config.verbose) {
      pw.println(progress + "\t" + cm)
      pw.println("## Typing complete")
    }
    pwFlanks.close

    val flanks = Source.fromFile(flanksOutput).getLines.toList.filterNot(_(0) == '#').map(_.drop(2))

    val cm2 = new CountMap[String]
    pw.println("##Aggregate typing with " + config.flankSize + "nt flanks")
    for (flank <- flanks.filter(_.size >= config.flankSize)) {
      cm2.count(flank.take(config.flankSize))

    }
    val values = cm2.map(_._2).toList
    if (config.verbose) {
      pw.println(cm2.mkString("\n"))

      pw.println(cm2.totalCount())

      pw.println("Min = " + values.min)
      pw.println("Max = " + values.max)
      val sum = values.foldLeft(0)(_ + _)
      pw.println("Avg. = " + (sum / values.size))
    }
    pw.println("##Filtered with " + config.threshold + "x threshold: ")

    pw.println(cm2.zipWithIndex.filter(p => p._1._2 >= config.threshold).map(zip => "> flank_" + zip._2 + ",support=" + zip._1._2+",S="+zip._1._1.take(1) + "\n" + zip._1._1.drop(1)).mkString("\n"))

    pw.println("# total supporting reads=" + cm2.filter(p => p._2 >= config.threshold).foldLeft(0)(_ + _._2))
    pw.println("# flanks found=" + cm2.filter(p => p._2 >= config.threshold).size)

    pw.close

  }

}