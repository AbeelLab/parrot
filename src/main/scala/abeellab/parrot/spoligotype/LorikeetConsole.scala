package abeellab.parrot.spoligotype

import abeellab.parrot.Main

object LorikeetConsole extends Main{

  override def main(args: Array[String]): Unit = {

    if (args.length == 0) {

      listInstructions
    } else {
      args(0) match {
        case "list" => listInstructions
        case "help" => listInstructions

        case "spoligotype" => LorikeetSpoligotyper.main(args.drop(1))
        
        case "multi-type" => MultiTyping.main(args.drop(1))
        
        case "_" => listInstructions
        case _ => listInstructions
      }
    }

  }

  def listInstructions() {
    println("Usage:java -jar lorikeet.jar [instruction] [instruction options...]")
    println("Instructions:")
    println("\tspoligotype            Spoligotype BAM file based on digital spoligotyping")
    println("\tmulti-type             Merge multiple spoligotype files together in a single file, renormalizing across multiple libraries when needed.")
    

  }

}