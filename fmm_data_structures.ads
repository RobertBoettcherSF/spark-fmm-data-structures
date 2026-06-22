-- FMM Data Structures Specification
-- Version: 0.08
-- Date: 2024-06-21
-- Description: SPARK formal specification for FMM data structures and algorithms

pragma SPARK_Mode (On);

package Fmm_Data_Structures is

   -- Core constants and types
   Max_Level      : constant := 8;
   Max_Boxes      : constant := 8 ** Max_Level;
   Max_Particles  : constant := 1_000_000;

   type Particle_Index is range 1 .. Max_Particles;
   type Box_Index is range 0 .. Max_Boxes - 1;
   type Level is range 0 .. Max_Level;
   type Morton_Index is mod 2**64;

   type Vector_3D is record
      X, Y, Z : Float;
   end record;

   type Particle_Array is array (Particle_Index) of Vector_3D;
   type Sort_Record is record
      Box_Id : Box_Index;
      Rank   : Natural;
   end record;
   type Sort_Array is array (Particle_Index) of Sort_Record;
   type Histogram_Array is array (Box_Index) of Natural;
   type Bookmark_Array is array (Box_Index) of Natural;
   type Non_Empty_Box_Array is array (Positive range <>) of Box_Index;

   -- Algorithm 1: Parallel-Pseudo-Sort (Fixed-Grid-Method)
   procedure Pseudo_Sort
     (Particles : in  Particle_Array;
      Sort_Idx  : out Sort_Array;
      Bin       : out Histogram_Array)
     with
       Global => null,
       Depends => (Sort_Idx => Particles, Bin => Particles),
       Pre    => (for all I in Particle_Index => Particles(I).X >= 0.0 and Particles(I).Y >= 0.0 and Particles(I).Z >= 0.0),
       Post   => (for all I in Box_Index => Bin (I) <= Max_Particles);

   -- Algorithm 3: GET-BOOKMARK-AND-BOX-INDEX
   procedure Compute_Bookmarks
     (Bin          : in  Histogram_Array;
      Bookmark     : out Bookmark_Array;
      Non_Empty_Id : out Non_Empty_Box_Array)
     with
       Global => null,
       Depends => (Bookmark => Bin, Non_Empty_Id => Bin);

   -- Morton index computation (3D bit interleaving)
   function Compute_Morton_Index (P : Vector_3D; L : Level) return Morton_Index
     with
       Pre  => P.X >= 0.0 and P.X <= Float(Max_Boxes) and
               P.Y >= 0.0 and P.Y <= Float(Max_Boxes) and
               P.Z >= 0.0 and P.Z <= Float(Max_Boxes);

end Fmm_Data_Structures;
