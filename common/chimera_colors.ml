(* Top 20 colors supported by Chimera:
   cf. http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/colortables.html

   Chimera's BILD format specification:
   http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/bild.html *)

type color =
  | Red
  | Orange_red
  | Orange
  | Yellow
  | Green
  | Forest_green
  | Cyan
  | Light_sea_green
  | Blue
  | Cornflower_blue
  | Medium_blue
  | Purple
  | Hot_pink
  | Magenta
  | White
  | Light_gray
  | Gray
  | Dark_gray
  | Dim_gray
  | Black

let to_bild = function
  | Red             -> ".color red"
  | Orange_red      -> ".color orange red"
  | Orange          -> ".color orange"
  | Yellow          -> ".color yellow"
  | Green           -> ".color green"
  | Forest_green    -> ".color forest green"
  | Cyan            -> ".color cyan"
  | Light_sea_green -> ".color light sea green"
  | Blue            -> ".color blue"
  | Cornflower_blue -> ".color cornflower blue"
  | Medium_blue     -> ".color medium blue"
  | Purple          -> ".color purple"
  | Hot_pink        -> ".color hot pink"
  | Magenta         -> ".color magenta"
  | White           -> ".color white"
  | Light_gray      -> ".color light gray"
  | Gray            -> ".color gray"
  | Dark_gray       -> ".color dark gray"
  | Dim_gray        -> ".color dim gray"
  | Black           -> ".color black"

(* useful to "rotate" colors *)
let next = function
  | Red             -> Orange_red
  | Orange_red      -> Orange
  | Orange          -> Yellow
  | Yellow          -> Green
  | Green           -> Forest_green
  | Forest_green    -> Cyan
  | Cyan            -> Light_sea_green
  | Light_sea_green -> Blue
  | Blue            -> Cornflower_blue
  | Cornflower_blue -> Medium_blue
  | Medium_blue     -> Purple
  | Purple          -> Hot_pink
  | Hot_pink        -> Magenta
  | Magenta         -> White
  | White           -> Light_gray
  | Light_gray      -> Gray
  | Gray            -> Dark_gray
  | Dark_gray       -> Dim_gray
  | Dim_gray        -> Black
  | Black           -> Red

(* some BILD commands:
   .sphere x y z r
   .comment ...
   .transparency 0.90
   .box x1 y1 z1 x2 y2 z2
*)
