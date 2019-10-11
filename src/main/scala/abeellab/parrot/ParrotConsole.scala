package abeellab.parrot

import atk.util.Tool
import abeellab.parrot.spoligotype.LorikeetConsole
import abeellab.parrot.rflp.RFLPTyper
import abeellab.parrot.spoligotype.MultiTyping
import abeellab.parrot.spoligotype.LorikeetSpoligotyper
import abeellab.macaw.MacawSNPtyper



trait Main extends Tool {
  def main(args: Array[String]) {}

}

object ParrotConsole extends Tool {

  def getDeclaredFields(cc: AnyRef) = {
    val m = (Map[String, Any]() /: cc.getClass.getDeclaredFields) { (a, f) =>
      f.setAccessible(true)
      a + (f.getName -> f.get(cc))
    }
    m.toList.sortBy(_._1)
  }

  val instructions: Map[String, Main] = Map(
    "rflp-is6110" -> RFLPTyper,
    "spoligotype" -> LorikeetSpoligotyper,
    "spoligo-multitype" -> MultiTyping,
    "macaw" ->MacawSNPtyper
    )
  def main(args: Array[String]): Unit = {

    if (args.length == 0 || !instructions.contains(args(0))) {

      listInstructions
    } else {
      val name = args(0)
      val obj: Main = instructions(name)
      obj.main(args.drop(1))

    }

  }

  def listInstructions() {
    println("Usage:java -jar parrot.jar <instruction> [options...]")
    println("Instructions:")
    println(instructions.toList.sortBy(_._1).map(f => String.format("    %1$-20s", f._1) + f._2.description).mkString("\n"))

  }

}