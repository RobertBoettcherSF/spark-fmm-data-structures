-- FMM Data Structures Implementation
-- Version: 0.08
-- Date: 2024-06-21
-- Description: SPARK implementation of FMM data structures and algorithms

pragma SPARK_Mode (On);
with Fmm_Data_Structures; use Fmm_Data_Structures;

package body Fmm_Data_Structures is

   -- Bit interleaving for 3D Morton index
   function Interleave_Bits (X_In, Y_In, Z_In : Integer) return Morton_Index is
      Result : Morton_Index := 0;
      Temp_X : Integer := X_In;
      Temp_Y : Integer := Y_In;
      Temp_Z : Integer := Z_In;
   begin
      for I in 0 .. 20 loop
         Result := Result or
           (Morton_Index (Temp_X mod 2) * 2**(3*I)) or
           (Morton_Index (Temp_Y mod 2) * 2**(3*I + 1)) or
           (Morton_Index (Temp_Z mod 2) * 2**(3*I + 2));
         Temp_X := Temp_X / 2;
         Temp_Y := Temp_Y / 2;
         Temp_Z := Temp_Z / 2;
      end loop;
      return Result;
   end Interleave_Bits;

   -- Morton index for a particle at a given level
   function Compute_Morton_Index (P : Vector_3D; L : Level) return Morton_Index is
      Scale : constant Float := Float (2 ** Integer(L));
      X_Quantized : Integer := Integer (P.X * Scale);
      Y_Quantized : Integer := Integer (P.Y * Scale);
      Z_Quantized : Integer := Integer (P.Z * Scale);
   begin
      return Interleave_Bits (X_Quantized, Y_Quantized, Z_Quantized);
   end Compute_Morton_Index;

   -- Algorithm 1: Parallel-Pseudo-Sort
   procedure Pseudo_Sort
     (Particles : in  Particle_Array;
      Sort_Idx  : out Sort_Array;
      Bin       : out Histogram_Array)
   is
   begin
      Bin := (others => 0);
      Sort_Idx := (others => (Box_Id => 0, Rank => 0));
      for I in Particle_Index loop
         declare
            Box_Id : Box_Index := Box_Index (Compute_Morton_Index (Particles (I), Max_Level) mod Max_Boxes);
         begin
            Sort_Idx (I).Box_Id := Box_Id;
            Sort_Idx (I).Rank   := Bin (Box_Id);
            Bin (Box_Id)        := Bin (Box_Id) + 1;
         end;
      end loop;
   end Pseudo_Sort;

   -- Algorithm 3: GET-BOOKMARK-AND-BOX-INDEX
   procedure Compute_Bookmarks
     (Bin          : in  Histogram_Array;
      Bookmark     : out Bookmark_Array;
      Non_Empty_Id : out Non_Empty_Box_Array)
   is
      Current_Sum : Natural := 0;
      Non_Empty_Count : Natural := 0;
   begin
      -- Count non-empty boxes and compute prefix sum
      for I in Box_Index loop
         pragma Loop_Invariant (Non_Empty_Count <= Natural(I) + 1);
         pragma Loop_Invariant (Current_Sum <= Max_Particles and Current_Sum + Bin(I) <= Max_Particles);
         if Bin (I) > 0 then
            Non_Empty_Count := Non_Empty_Count + 1;
         end if;
         Bookmark (I) := Current_Sum;
         Current_Sum := Current_Sum + Bin (I);
      end loop;

      -- Populate Non_Empty_Id
      Non_Empty_Id := (1 .. Non_Empty_Count => 0);
      declare
         J : Positive := 1;
      begin
         for I in Box_Index loop
            pragma Loop_Invariant (J <= Non_Empty_Count + 1 and J >= 1);
            if Bin (I) > 0 then
               Non_Empty_Id (J) := I;
               J := J + 1;
            end if;
         end loop;
      end;
   end Compute_Bookmarks;

end Fmm_Data_Structures;
