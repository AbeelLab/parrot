package abeellab.macaw

import net.sf.samtools.SAMRecordIterator
import net.sf.samtools.BAMFileReader
import scala.collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import java.io.File
import be.abeel.util.CountMap
import java.io.PrintWriter
import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult
import net.sf.jannot.utils.SequenceTools
import net.sf.jannot.refseq.Sequence
import net.sf.jannot.refseq.MemorySequence
import scala.io.Source
import atk.util.Tool

import java.text.NumberFormat
import java.util.Locale

import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils
import abeellab.parrot.Main
import atk.io.NixWriter
import edu.northwestern.at.utils.math.statistics.FishersExactTest
import java.util.HashMap

/**
 *
 * SNP typer
 * Features:
 * - 100% matches only
 * -runs entire BAM file
 *
 */
object MacawSNPtyper extends Main {

  override def description = "Determines the lineage(s) of a MTB sequencing sample. Reports relative abundance of multiple lineages if present, i.e. mixed infections. Typing is done through a marker set."
  
  override val version = """
    2015/01/16:    Initial release
    2015/03/03:    Changed output logic to only output a single marker type for all files by default
		  		   Added option to revert to the old behavior
    2019/09/11:    Clean up console output, integrated post-processing in primary tool. Merged with other typing tools.
    """

  def revcomp(read: Array[Byte]) = {
    val out = Array.ofDim[Byte](read.length)
    for (i <- 0 until read.length) {
      out(read.length - i - 1) = SequenceTools.complement(read(i).toChar).toByte
    }
    out
  }
  case class Config(val detailed: Boolean = false, val markerFile: String = null, val outputFile: String = null, files: List[File] = List(), val threshold: Int = 5)
  /**
   * args(0) = output file
   *
   *
   */
  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar macaw.jar") {
      opt[String]("marker") action { (x, c) => c.copy(markerFile = x) } text ("File containing marker sequences. This file has to be a multi-fasta file with the headers indicating the name of the markers.") //, { v: String => config.spacerFile = v })
      opt[String]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]('t', "threshold") action { (x, c) => c.copy(threshold = x) } text ("Threshold to determine absence or presence of a marker (default=5)")
      //opt[Unit]("detailed") action { (_, c) => c.copy(detailed = true) } text ("Output digital marker types per input file. (default=false) ")
      arg[File]("<file>...") unbounded () required () action { (x, c) => c.copy(files = c.files :+ x) } text ("input files")

    }
    parser.parse(args, Config()) map { config =>
      /* Load spacers */
      val lines = if (config.markerFile != null) tLines(config.markerFile).toList else scala.io.Source.fromInputStream(MacawSNPtyper.getClass().getResourceAsStream("/subset3LongMarkers.txt")).getLines().filterNot(f => f.startsWith("#") || f.trim.size == 0).toList;
//      assume(config.outputFile != null)
      val of=if(config.outputFile ==null) new File(config.files(0)+".macaw") else  config.outputFile
      val pw = new NixWriter(of+".intermediate")

      pw.println(generatorInfo)

      pw.println("# Marker input: " + (if (config.markerFile == null) "(default)" else config.markerFile))

      val in = (lines.grouped(2).map(f => (f(0).substring(1), f(1))).toList)

      val repeatSequences = List(("repeat", "GTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGA"), ("left_repeat", "GTTTCCGTCCCC"), ("middle_repeat", "TCTCGGGGTTTT"), ("right_repeat", "GGGTCTGACGA"))

      val forwardSpacers = in ++ repeatSequences // ++ mirus()

      val rcSpacers = forwardSpacers.map { case (x, y) => ("RC_" + x, new String(revcomp(y.getBytes))) }

      val spacers = forwardSpacers ++ rcSpacers

      pw.println("# Input files: " + config.files)

      /* Prep index */
      val tree = new AhoCorasick();
      for ((id, spacer) <- spacers) {
        tree.add(spacer.getBytes(), id)
      }
      tree.prepare();

      val nf = NumberFormat.getInstance(Locale.US)
      nf.setMaximumFractionDigits(6)
      var cm = new CountMap[String]
      for (inputFile <- config.files) {
        pw.println("# Processing: " + inputFile)
        /* Connect to bam file*/
        val inputSam = new SAMFileReader(inputFile);

        if (config.detailed)
          cm = new CountMap[String]

        /* Iterate over bamfile */
        val it: SAMRecordIterator = inputSam.iterator()
        var progress: Int = 0
        var totalCoverage: Long = 0

        val time = System.currentTimeMillis();
        val nano = System.nanoTime()

        while (it.hasNext()) {
          val sr = it.next()

          val read = sr.getReadBases()

          totalCoverage += sr.getReadLength()

          val result = tree.search(read);

          for (a <- result) {
            for (s <- a.asInstanceOf[SearchResult].getOutputs()) {
              cm.count(s.asInstanceOf[String])
            }
          }

          progress += 1
          if (progress % 100000 == 0) {
            print(".")

          }
        }
        println

        /* close bam file */
        it.close

        //        if (config.detailed)
        //          output(pw, spacers, cm, config)

      }
      if (!config.detailed)
        output(pw, spacers, cm, config)
      pw.close

      /**
       * Interpret marker file
       */
      val linCountsPresent = new CountMap[String] //number present in the indicated group
      val linCountsPresentInverse = new CountMap[String] //# number present in all except the indicated group
      val linCountsAbsent = new CountMap[String] //		# number absent in the indicated group
      val linCountsAbsentInverse = new CountMap[String] //	# number absent in all except the indicated group
      val totalMarkersPerLin = new CountMap[String] //		# total number of markers in each lineage
      val markerCounts = new CountMap[String] //	# read depths for markers in each lineage
      val avgMarkerPerLin = new HashMap[String,Double] //	# average read depth for markers in each lineage
      val linPercentTotal = new HashMap[String,Double]
      val linPercentEnd = new CountMap[String]
      
      //predLin = []					# predicted lineages

      val fileContent=tLines(of+".intermediate")
      for (line <- fileContent) {
        val arr=line.split("\t")
        val mCount=arr(1).toInt
        val lin=arr(0).split("_")(2)
        totalMarkersPerLin.count(lin)
        markerCounts.count(lin,mCount)
        arr(3) match{
          case "present" =>linCountsPresent.count(lin)
          case "absent" =>linCountsAbsent.count(lin)
        }
        
        
        

      }
      val sumAllPresent=linCountsPresent.values().toList.map(_.toInt).sum
      val sumAllAbsent=linCountsAbsent.values().toList.map(_.toInt).sum
      
      val allKeys=(linCountsPresent.keySet()++linCountsAbsent.keySet()).toList.sorted
      /* Calculate p-values */
      val fisherPValues = new HashMap[String,Double] //		# p values for each lineage
      for(key<-allKeys){
          linCountsPresentInverse.count(key,sumAllPresent - linCountsPresent(key))
          linCountsAbsentInverse.count(key,sumAllAbsent - linCountsAbsent(key))
//          println(linCountsPresent(key),  linCountsAbsent(key),linCountsPresentInverse(key), linCountsAbsentInverse(key))
         val rawP=FishersExactTest.fishersExactTest(linCountsPresent(key),  linCountsPresentInverse(key),linCountsAbsent(key), linCountsAbsentInverse(key))
         
//         println(rawP.toList)
         val correctedP= rawP(2) * allKeys.size
		     fisherPValues(key) = if(correctedP>1)1 else correctedP
		     avgMarkerPerLin(key) = 1.0*markerCounts(key) / totalMarkersPerLin(key)

      }
      val allLinDepth=avgMarkerPerLin.values().sum
      for(key<-allKeys){
        linPercentTotal(key) = avgMarkerPerLin(key) / allLinDepth
      }
      
      
      val predicted=allKeys.filter(fisherPValues(_)<=0.05)
      
      val pw2=new NixWriter(of.toString)
      predicted.size match{
        case 0=>
		      pw2.println("Lineages:\tNONE")
		      pw2.println("Mixed:\tFALSE")  
		    case 1=>
		      pw2.println("Lineages:\t"+predicted.head)
		      pw2.println("Mixed:\tFALSE")
        case _=>
          val totalLinDepth=predicted.map(f=>avgMarkerPerLin(f)).sum
		      val output=for(lineage<-predicted)yield{
		        val linPercent = avgMarkerPerLin(lineage) / totalLinDepth
		        lineage+" ("+linPercent+")"
		      }
          pw2.println("Lineages:\t"+output.mkString(","))
          pw2.println("Mixed:\tTRUE")
      }
      
		
      pw2.println("##lineage\tpresent markers\tabsent markers\traw markers\tfraction\tcorrected p-value")
      for(key<-allKeys){
         
         pw2.println(key+"\t"+linCountsPresent(key)+"\t"+linCountsAbsent(key)+"\t"+markerCounts(key)+"\t"+nf.format(linPercentTotal(key))+"\t"+fisherPValues(key))
      }
      pw2.close
  
      
      
    } getOrElse {
      println("Could not interpret command-line arguments, quitting!")
      System.exit(-1)
    }
  }

  def output(pw: PrintWriter, spacers: List[(String, String)], cm: CountMap[String], config: Config) {
    val listx = spacers.filter(p => !p._1.contains("repeat")).map(f => cm.get(f._1).toInt)

    pw.println("# number of spacers = " + spacers.size)

    val groupedSpacers = spacers.groupBy(pair => pair._1.replaceAll("RC_", ""))

    //    println("GS: " + groupedSpacers.mkString("\n"))
    val buffer = new StringBuffer()
    pw.println("# Marker\tread-depth\tp-value\tA/P")
    var idx = 0

//    println("KS: " + cm.keySet())
    for (gs <- groupedSpacers.filterNot(_._1.contains("repeat")).toList.sortBy(_._1)) {
      idx += 1
      assert(gs._2.size == 2)
//      println("GGS: " + gs)
      val z1 = (gs._2.map(p => cm.get(p._1).toInt)).toList
      val z = z1.sum

      pw.println(gs._1 + "\t" + z + "\t" + nf.format(if (z >= config.threshold) 0 else 1) + "\t" + (if (z >= config.threshold) "present" else "absent"))
      buffer.append(if (z >= config.threshold) "1" else "0")
      //      if (buffer.length() % 11 == 10)
      //        buffer.append(" ")

    }
    pw.println("## Digital markertype: \n#" + buffer.toString().grouped(10).mkString(" "))
    
  }
}
