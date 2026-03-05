/*
 *  This program read in the MODIS Specification (CDL format) files one line
 *  at a time.  It then parses the line according to the data type and
 *  creates two fortran subroutines get_spec_info.f and get_spec_data.f
 *  from the Specification information.  
 *
 *  Created by Walter Wolf 05/19/97
 *  University of Wisconsin - Madison  Space Science and Engineering Center
 *  
 *  If any problems, questions, or comments about this code, contact 
 *  Walter Wolf at the following email address:
 *
 *  wolf@frogman.ssec.wisc.edu
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct _StringHeaderRec {
   char  Strings[50][320];
   char  LastChar;
   int   NumOfStrings;
   int   EqualSign;
   int   SemiColon;
}  StringHeaderRec;
typedef struct _StringHeaderRec *StringHeader;

typedef struct _NewLineHeaderRec {
   char  Strings[6][320];
   int   Lines;
}  NewLineHeaderRec;
typedef struct _NewLineHeaderRec *NewLineHeader;


/*  Read in one line of data */

int getline(FILE *file, char *s, int lim)
{
   int c, i;

   for (i=0; i<lim-1 && (c=fgetc(file)) != EOF && c != '\n'; ++i)
      s[i]=c;
   if (c == '\n') {
      s[i]=c;
      ++i;
   }
   s[i]='\0';
   return i;
}

/*
 *  Take a line from the spec file and separate it according to the data
 *  type.
 */

void separate (char *buff, StringHeader STR) 
{
   int i;
   int length;
   int StartPoint;
   int Flag = 0;

  /*  Clear the array */

   for (i = 0; i < 50; i++)
      strncpy((char *) STR->Strings[i], " ", 319);

  /*  Initialize the variables */

   length = strlen(buff);
   STR->NumOfStrings = 0;
   STR->EqualSign = 0;
   STR->SemiColon = 0;

  /*  Check the characters in the string, one at a time */

   for (i = 0; i < length; i++) {

     /*  If it is a semi-colon or a return sequence, then break */

      if (buff[i] == ';' || buff[i] == '\n') {
         Flag = 0;
         if (buff[i] == ';') STR->SemiColon = 1;
         break;
      }

    /*  If it is an equal sign not in a string, then set the EqualSign 
        variable */

      else if (buff[i] == '=' && Flag != 3) {
         STR->EqualSign = 1;
         continue;
      }

     /*  If it is the first quote, then set flag to 2 */

      else if (buff[i] == '"' && Flag != 3)
         Flag = 2;

    /*  Save the last character is it is not whitespace  */

     if (buff[i] != ' ' && buff[i] != '\t')
        STR->LastChar = buff[i];

     /*  if it is any one of the following characters, then break  */

/*
      if (i == 0 && (buff[i] == '\\' || buff[i] == '}' || buff[i] == '/'))
*/
      if (i == 0 && (buff[i] == '\\' || buff[i] == '}'))
         break;

      if (!Flag) {

        /*  Check white space or tab characters */

         if (buff[i] != ' ' && buff[i] != '\t') {
            StartPoint = i;
            Flag = 1;
         }
      }
      if (Flag == 1) {

        /*  While we are reading a string, then check for white space.
            If it exists, then the string is done */

         if (buff[i] == ' ' || buff[i] == '\t') {
            Flag = 0; 
            STR->NumOfStrings += 1;
         }
         else
            STR->Strings[STR->NumOfStrings][i-StartPoint] = buff[i];
      }
      else if (Flag == 2) {

        /*  Used to set up the character string read  */

         StartPoint = i + 1;
         Flag = 3;
      }
      else if (Flag == 3) {

        /*  Check for the end quote while reading in a string.  If the
            character is the second quote, then continue, else store the
            character */

         if (buff[i] == '"') {
            STR->NumOfStrings += 1;
            Flag = 0;
         }
         else {
            STR->Strings[STR->NumOfStrings][i-StartPoint] = buff[i];
         }
      }
   }
}

/*  
 *  trim_line takes a max 320 character line and breaks it up into
 *  fortran happy shorter strings
 */

void trim_line (StringHeader STR, NewLineHeader NLH)
{
   int i;
   int length;
   int lastchar, newchar=0;
   int startchar;

   for (i = 0; i < 6; i++)
      strncpy((char *) NLH->Strings[i], " ", 319);

   length = strlen(STR->Strings[0]);
   lastchar = 0;

  /*  Find the last non blank or non (\n) character */

   for (i = 0; i < length; i++) {
      if (STR->Strings[0][i] != ' ' && STR->Strings[0][i] != '\n' &&
          STR->Strings[0][i] != '\t') {
         if (STR->Strings[0][i] == '\\' && STR->Strings[0][i+1] == 'n')
            break;
        lastchar = i;
      }
   }

  /*  
   * This number 70 assumes that we are using the f77 option -col120
   *  (120 columns) in the compilation of the get_spec_info.f routine.
   */

   if (lastchar > 70) {
      
     /*  break the line in half */

     /* look for a clean break, a : or a - */

      newchar = 0;
      for (i = 0; i <= lastchar; i++) {
         if (STR->Strings[0][i] == ':' || STR->Strings[0][i] == '-')
            newchar = i;
      }

     /*  if no clean break, then find a blank space after the 30th character */

      if (lastchar > 140) {
         startchar = 80;
      }
      else {
         startchar = 60;
      }
      if (newchar > startchar || newchar < 15) {
         for (i = startchar; i <= lastchar; i++) {
            if (STR->Strings[0][i] == ' ') {
               newchar = i;
               break;
            }
         }
      }
      if (newchar == 0) newchar = startchar / 2;
   }

  /*  Copy the string(s) to the output strings */

   if (newchar > 0 && newchar < 1100) {
      for (i = 0; i <= newchar; i++) {
         NLH->Strings[0][i] = STR->Strings[0][i];
         NLH->Lines = 1;
      }
      for (i = newchar+1; i <= lastchar; i++) {
         NLH->Strings[1][i-newchar] = STR->Strings[0][i];
         NLH->Lines = 2;
      }
   }
   else {
      for (i = 0; i <= lastchar; i++) {
         NLH->Strings[0][i] = STR->Strings[0][i];
         NLH->Lines = 1;
      }
   }
}

/*
 *  Take one line read in and determine if there is anything but
 *  white space in this line
 */

int check_whitespace (char *buff)
{
   int i;
   int length, Start;

   length = strlen(buff);
   Start = -1;

   for (i = 0; i < length; i++) {
      if (buff[i] != ' ' && buff[i] != '\t' && buff[i] != '\n') {
         Start = i;
         break;
      }
   }
   return (Start);
}

/*
 *   Write a certain number of blanks to the output file
 *   This is used to set up the fortran code.
 */

void write_blanks (FILE *Out, int Number)
{
   int i;
  
   for (i = 0; i < Number; i++) 
      fprintf (Out, " ");
}


/*  Main program */

int main (int argc, char *argv[])
{
   FILE *InputFile, *OutputFile, *TextFile;
   char  FileName[320];
   char  buffer[320];
   char  hold[320];
   char  DimList[320];
   char  DimListVariables[5][320];
   char  MODName[15]; /*Increased by G.Britzolakis on 12/27/11 to 15 characters from 10 to fix bug*/
   char  SpecName[500][320];
   static char  SpecVar[1000][20][320];
   static char  Dimensions[1000][320];
   static char  Variables[1000][320];
   static char  OutLines[10000][320];
   static char  OutVariables[10000][320];
   static char  Prog2Lines[20000][320];
   int   OutSpaces[10000];
   int   Prog2Count = 0;
   int   Flag;
   int   i, j;
   int   DimCount, VarCount, SpecCount = 0;
   int   Dim2Count = 0, DataCount, DataVariables = 0;
   int   EqualFlag[10000];
   int   LineNumber, NumberOfLines;
   int   VariableMatch = 0, NumberOfVariables;
   int   Length, Spaces = 7;
   int   TitleCount = 0;
   int   HistoryCount = 0;
   int   DescCount = 0;
   int   DimFlag = 0;
   int   DimListCount = 0;
   int   DimListVarCount = 0;
   int   Global_String = 0;
   int   Global_Length = 0;
   long  Value = 0;
   static StringHeaderRec STR;
   static NewLineHeaderRec NLH;

   if (argc > 2) {
      printf("Too many input arguments!!\n");
      exit(0);
   }

  /*  Open up the output file */

   if ((OutputFile = fopen("get_spec_info.f","w")) == NULL) {
      printf ("Could not open the output file!!\n");
      exit(0);
   }

  /*  Write the top of the fortran subroutine */

   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "SUBROUTINE GET_SPEC_INFO (PROGRAMNAME, INPUTNAME,\n");
   Spaces += 17;
   write_blanks (OutputFile, 5);
   fprintf (OutputFile, "$ ");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "LONG_NAME, DIMLIST, UNITS_STR, VALID_RANGE, FILL_VALUE,\n");
   write_blanks (OutputFile, 5);
   fprintf (OutputFile, "$ ");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "NUMBERTYPE, SCALE_FACTOR, ADD_OFFSET, PARAM_TYPE,\n");
   write_blanks (OutputFile, 5);
   fprintf (OutputFile, "$ ");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "CELL_ACROSS_SWATH_SAMPLING, CELL_ALONG_SWATH_SAMPLING,\n");
   write_blanks (OutputFile, 5);
   fprintf (OutputFile, "$ ");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "DESCRIPTION, DESCCOUNT, TITLE, TITLECOUNT,\n");
   write_blanks (OutputFile, 5);
   fprintf (OutputFile, "$ ");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "HISTORY, HISTORYCOUNT, GEOLOCATION_POINTER)\n");

   Spaces -= 17;
   fprintf (OutputFile, "\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "implicit none\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "save\n\n");


   fprintf (OutputFile, "c---------------------------------------------------------------------\n");
   fprintf (OutputFile, "C!F77\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!DESCRIPTION:\n");
   fprintf (OutputFile, "c   Subroutine which extracts the correct file spec information\n");
   fprintf (OutputFile, "c   and creates an HDF output file for MOD35, MOD07, or MOD06.\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!INPUT PARAMETERS:\n");
   fprintf (OutputFile, "c   ProgramName   Ouput hdf file name\n");
   fprintf (OutputFile, "c   InputName     Name of Input SDS\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!OUTPUT PARAMETERS:\n");
   fprintf (OutputFile, "c   long_name                   HDF long name of SDS\n");
   fprintf (OutputFile, "c   dimlist                     Dimension list for SDS\n");
   fprintf (OutputFile, "c   units_str                   Unit of current SDS\n");
   fprintf (OutputFile, "c   valid_range                 Valid range of current SDS\n");
   fprintf (OutputFile, "c   Fill_Value                  Missing data value for given SDS\n");
   fprintf (OutputFile, "c   NumberType                  Type of data for current SDS\n");
   fprintf (OutputFile, "c   scale_factor                Scale factor for current SDS\n");
   fprintf (OutputFile, "c   add_offset                  Add/offset values for current SDS\n");
   fprintf (OutputFile, "c   Param_Type                  Parameter type of current SDS\n");
   fprintf (OutputFile, "c   Cell_Across_Swath_Sampling  Sampling for current SDS\n");
   fprintf (OutputFile, "c   Cell_Along_Swath_Sampling   Sampling for current SDS\n");
   fprintf (OutputFile, "c   description                 SDS description\n");
   fprintf (OutputFile, "c   DescCount                   Number of lines to be included in description\n");
   fprintf (OutputFile, "c   title                       SDS title\n");
   fprintf (OutputFile, "c   TitleCount                  Number of title lines\n");
   fprintf (OutputFile, "c   history                     SDS history summary information\n");
   fprintf (OutputFile, "c   HistoryCount                Number of history lines\n");
   fprintf (OutputFile, "c   Geolocation_Pointer         SDS geolocation pointer\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!REVISION HISTORY:\n");
   fprintf (OutputFile, "c   This FORTRAN source file is created by spec2fort.c\n");
   fprintf (OutputFile, "c   spec2fort.c was designed by Walter.Wolf@ssec.wisc.edu\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!TEAM-UNIQUE HEADER:\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!REFERENCES AND CREDITS:\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!END\n");
   fprintf (OutputFile, "c---------------------------------------------------------------------\n\n");

   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "include \'mapi.inc\'\n\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "include \'hdf.inc\'\n\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "include \'dffunc.inc\'\n\n");

  /*  Declare the fortran input variable types  */

   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "character*(*) ProgramName, InputName, Geolocation_Pointer\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "character*(*) long_name, units_str, Param_Type, dimlist\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "double precision valid_range(*), Fill_Value\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "double precision scale_factor, add_offset\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "integer Cell_Across_Swath_Sampling(*)\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "integer Cell_Along_Swath_Sampling(*)\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "character*(*) title, history\n");
   write_blanks (OutputFile, Spaces);
/*
   fprintf (OutputFile, "character*(*) description\n");
*/
   fprintf (OutputFile, "character*(*) description(*)\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "integer TitleCount, HistoryCount, DescCount, NumberType\n\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "integer unlimited\n");

  /*  
   *  Open up the input file.  This input file contains a list of
   *  CDL files with names like  MOD05.V2.CDL
   */

   if ((TextFile = fopen(argv[1],"r")) == NULL) {
      printf ("Could not open the input file!!\n");
      exit(0);
   }

  /*  Read in the spec file name */

   NumberOfLines = 0;
   NumberOfVariables = 0;
   while (fscanf (TextFile, "%s", (char *) &FileName) > 0) {
   
     /*  Open up each CDL spec file one at a time  */

      if ((InputFile = fopen(FileName,"r")) == NULL) {
         printf ("Could not open the input file!!\n");
         exit(0);
      }

     /*  write the if statement in the fortran code for the spec info for
         each MODIS program name */

     /* 
      * Determine the number of characters in the name before 
      * the first period
      */

      if (NumberOfLines == 0) {
         OutSpaces[NumberOfLines] = Spaces;
         sprintf (OutLines[NumberOfLines], "unlimited = 0\n");
         sprintf (Prog2Lines[Prog2Count], "unlimited = 0\n");
         Prog2Count += 1;
         NumberOfLines += 1;
         OutSpaces[NumberOfLines] = Spaces;
         sprintf (OutLines[NumberOfLines], "TitleCount = 0\n");
         NumberOfLines += 1;
         OutSpaces[NumberOfLines] = Spaces;
         sprintf (OutLines[NumberOfLines], "HistoryCount = 0\n");
         NumberOfLines += 1;
      }

      i = strcspn((char *) FileName, ".");

      strncpy (MODName, FileName, i);
      OutSpaces[NumberOfLines] = Spaces;
      sprintf (OutLines[NumberOfLines], "if (ProgramName .EQ. \"%s\") then\n",
               MODName);
      sprintf (Prog2Lines[Prog2Count], "if (ProgramName .EQ. \"%s\") then",
               MODName);
      Prog2Count += 1;
      NumberOfLines += 1;
      Spaces += 3;

     /*  Clear out a character array  */

      for (i = 0; i < 50; i++) 
         strncpy((char *) SpecName[SpecCount], " ", 319);

     /*  Initialize variables */

      Flag = 0;
      DimCount = 0;
      VarCount = 0;
      SpecCount = 0;
      DescCount = 0;
      LineNumber = -1;

   /*  read in the data one line at a time and then start parsing */

      while (getline(InputFile, buffer, 320) > 0) {
         LineNumber += 1;
         EqualFlag[LineNumber] = 0;
         if (check_whitespace(buffer) < 0) {
            Flag = 0;
            continue;
         }
 
         switch (Flag) {

         case 0:

           /*  Write the total number of history, title, and/or 
               description lines if they exist */

            if (TitleCount > 0) {
/*
               sprintf (OutLines[NumberOfLines],
                        "TitleCount = %d \n", TitleCount);
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
*/
               TitleCount = 0;
            }
            if (HistoryCount > 0) {
/*
               sprintf (OutLines[NumberOfLines],
                        "HistoryCount = %d \n", HistoryCount);
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
*/
               HistoryCount = 0;
            }
            if (Dim2Count > 0) {
               sprintf (Prog2Lines[Prog2Count], "   DimCount = %d ", Dim2Count);
               Prog2Count += 1;
               Dim2Count = 0;
            }
            if (DataVariables > 0) {
               sprintf (Prog2Lines[Prog2Count], "   NumOfData = %d ", 
                        DataVariables);
               Prog2Count += 1;
               DataVariables = 0;
            }
            if (DescCount > 0) {
               sprintf (OutLines[NumberOfLines],
                        "DescCount = %d \n", DescCount);
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
               DescCount = 0;
            }

           /*  Read in the next line and determine what type of fortran
               data to write out */

            separate(buffer, &STR);
            if (strncmp((char *) STR.Strings[0], "dimensions:", 11) == 0) {
               Flag = 1;
               Dim2Count = 0;
             /*
              *  These two lines are commented out since they are not needed
              *  in the get_spec_info subroutine.
              *
               sprintf (OutLines[NumberOfLines], "\n");
               NumberOfLines += 1;
              */
               OutSpaces[NumberOfLines] = Spaces;
               getline(InputFile, buffer, 320);
               LineNumber += 1;
            }
            else if (strncmp((char *) STR.Strings[0], "variables:", 10) == 0) {
               Flag = 2;
               sprintf (OutLines[NumberOfLines], "\n");
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
               getline(InputFile, buffer, 320);
               LineNumber += 1;
            }
/*
            else if (strncmp((char *) STR.Strings[1], ":title", 6) == 0) {
               Flag = 4;
               sprintf (OutLines[NumberOfLines], "\n");
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
               getline(InputFile, buffer, 320);
               LineNumber += 1;
            }
            else if (strncmp((char *) STR.Strings[1], ":history", 8) == 0) {
               Flag = 5;
               sprintf (OutLines[NumberOfLines], "\n");
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
               getline(InputFile, buffer, 320);
               LineNumber += 1;
            }
*/
            else if (strncmp((char *) STR.Strings[0], "data:", 5) == 0) {
               Flag = 6;
               getline(InputFile, buffer, 320);
               LineNumber += 1;
            }
            else if (strcmp((char *) STR.Strings[0], "byte")   == 0 ||
                     strcmp((char *) STR.Strings[0], "short")  == 0 ||
                     strcmp((char *) STR.Strings[0], "long")   == 0 ||
                     strcmp((char *) STR.Strings[0], "float")  == 0 ||
                     strcmp((char *) STR.Strings[0], "double") == 0) {
             
              /*  Close off the if statement if not the first time through */

               if (SpecCount > 0) {
                  Spaces -= 3;
                  OutSpaces[NumberOfLines] = Spaces;
                  sprintf (OutLines[NumberOfLines], "End If\n\n");
                  NumberOfLines += 1;
               }
               else {
                  OutSpaces[NumberOfLines] = Spaces;
                  sprintf (OutLines[NumberOfLines], "\n");
                  NumberOfLines += 1;
               }

               strncpy((char *) SpecName[SpecCount], " ", 319);

              /*  Open up the fortran if statement for the data types */

               if (STR.NumOfStrings >= 4) {
                  for (i = 0; i < STR.NumOfStrings; i++) {
                     if (STR.Strings[i][0] == '(') {
                        sprintf (OutLines[NumberOfLines],
                                "if (InputName .EQ. \"%s\") then\n",
                                 SpecName[SpecCount]);
                        OutSpaces[NumberOfLines] = Spaces;
                        Spaces += 3;
                        NumberOfLines += 1;
                        break;
                     }
                     else {
                        strcpy((char *) SpecName[SpecCount], STR.Strings[1]);
                     }
                  }
               }
               else {
                  for (i = 0; i < (int) strlen(STR.Strings[1]); i++) {
                     if (STR.Strings[1][i] == '(') {
                        sprintf (OutLines[NumberOfLines], 
                                "if (InputName .EQ. \"%s\") then\n",
                                 SpecName[SpecCount]);
                        OutSpaces[NumberOfLines] = Spaces;
                        Spaces += 3;
                        NumberOfLines += 1;
                        break;
                     }
                     else {
                        SpecName[SpecCount][i] = STR.Strings[1][i];
                     }
                  }
               }
               SpecCount += 1;
               Flag = 3;

              /*  Write the data type to the fortran file */

               OutSpaces[NumberOfLines] = Spaces;
               if (strcmp((char *) STR.Strings[0], "byte")   == 0)
                  sprintf (OutLines[NumberOfLines], "NumberType = DFNT_INT8\n");
               else if (strcmp((char *) STR.Strings[0], "short")   == 0)
                  sprintf (OutLines[NumberOfLines], 
                           "NumberType = DFNT_INT16\n");
               else if (strcmp((char *) STR.Strings[0], "long")   == 0)
                  sprintf (OutLines[NumberOfLines], 
                           "NumberType = DFNT_INT32\n");
               else if (strcmp((char *) STR.Strings[0], "float")   == 0)
                  sprintf (OutLines[NumberOfLines], 
                           "NumberType = DFNT_FLOAT32\n");
               else if (strcmp((char *) STR.Strings[0], "double")   == 0)
                  sprintf (OutLines[NumberOfLines], 
                           "NumberType = DFNT_FLOAT64\n");

               NumberOfLines += 1;

              /*  Write out the dimension list for the particular data type */

               strncpy((char *) DimList, " ", 319);
               for (i = 0; i < 5; i++) 
                  strncpy((char *) DimListVariables[i], " ", 319);

               strncpy((char *) DimList, " ", 319);
               if (STR.NumOfStrings >= 4) {
                  DimFlag = 0;
                  DimListCount = 0;
                  DimListVarCount = -1;
                  for (i = 0; i < STR.NumOfStrings; i++) {
                     if (STR.Strings[i][0] == '(' ) {
                        DimFlag = 1;
                        continue;
                     }
                     if (STR.Strings[i][0] == ')' )
                        break;
                     if (STR.Strings[i][0] == ',' ) {
                        DimFlag = 1;
                        continue;
                     }
                     if (DimFlag) {
                        DimListVarCount += 1;
                        strcpy (DimListVariables[DimListVarCount], 
                                STR.Strings[i]);
                        DimFlag = 0;
                     }
                  }
               }
               else {
                  Length = strlen(STR.Strings[1]);
                  DimFlag = 0;
                  DimListCount = 0;
                  DimListVarCount = 0;
                  for (i = 0; i < Length; i++) {
                     if (STR.Strings[1][i] == '(' ) {
                        DimFlag = 1;
                        continue;
                     }
                     if (STR.Strings[1][i] == ')' )
                        break;
                     if (STR.Strings[1][i] == ':' )
                        DimFlag = 0;
                     if (STR.Strings[1][i] == ',' ) {
                        DimFlag = 1;
                        DimListVarCount += 1;
                        DimListCount = 0;
                        continue;
                     }
                     if (DimFlag) {
                        DimListVariables[DimListVarCount][DimListCount] = 
                                         STR.Strings[1][i];
                        DimListCount += 1;
                     }
                  }
               }
            
              /*  In the CDL file, the dimension name for each variable is
               *  written in C order.  We have to write out the dimension
               *  names in Fortran order since we are creating a Fortran 
               *  subroutine.
               */

               OutSpaces[NumberOfLines] = Spaces;
               sprintf (DimList,
                      "DimList = \n     $ \"%s", 
                      DimListVariables[DimListVarCount]);
               for (i = DimListVarCount - 1; i >= 0; i--) {
                  if (strlen(DimListVariables[i]) > 40)
                     sprintf (DimList, "%s,\" //\n     $ \"%s", 
                              DimList, DimListVariables[i]);
                  else
                     sprintf (DimList, "%s,%s", DimList, DimListVariables[i]);
               }

               sprintf (DimList, "%s\"", DimList);
               sprintf (OutLines[NumberOfLines], "%s\n", DimList);
               NumberOfLines += 1;
            }
            break;

         case 1:

            Dim2Count += 1;

           /*  Read in the dimensioned variables */

            separate(buffer, &STR);

           /*  Store the Dimension definitions to an array */

           /*  This loop is to strip out the colon  after the
            *  Dimension variable Name */

            strncpy (Dimensions[DimCount], " ", 319);; 
            for (j = 0; j < (int) strlen(STR.Strings[0]); j++) {
               if (STR.Strings[0][j] == ':')
                  break;
               else
                  Dimensions[DimCount][j] = STR.Strings[0][j]; 
            }
            DimCount += 1;
      
           /*  Now copy the Number */

            strcpy(Dimensions[DimCount], STR.Strings[1]); 
            DimCount += 1;

           /*  Now write it out */

/*
 *  This section is commented out since these variables are not used in the
 *  subroutine.  They are used in the subroutine get_spec_data so they are
 *  not needed  here.
 *
            if (STR.EqualSign) {
               OutSpaces[NumberOfLines] = Spaces;
               sprintf (OutLines[NumberOfLines], "%s = %s\n", 
                        Dimensions[DimCount-2], 
                        Dimensions[DimCount-1]);
               NumberOfLines += 1;
 
               VariableMatch = 1;
               for (i = 0; i < NumberOfVariables; i++) {
                  if (strcmp (Dimensions[DimCount-2], OutVariables[i]) == 0) {
                     VariableMatch = 0;
                     break;
                  }
               }
               if (VariableMatch) {
                  write_blanks (OutputFile, 7);
                  fprintf (OutputFile, "integer %s\n", Dimensions[DimCount-2]);
                  strcpy (OutVariables[NumberOfVariables], 
                          Dimensions[DimCount-2]);
                  NumberOfVariables += 1;
               }
            }
*/
  
           /*  Save data for the second program */

            sprintf (Prog2Lines[Prog2Count], "   DimNames(%d) = \"%s\"", 
                     Dim2Count, Dimensions[DimCount-2]);
            Prog2Count += 1;
            sprintf (Prog2Lines[Prog2Count], "   DimNumbers(%d) = %s", 
                     Dim2Count, Dimensions[DimCount-1]);
            Prog2Count += 1;
            break;
         
         case 2:

           /*  Read in the variable definitions */

            separate(buffer, &STR);

            Global_Length = strlen(buffer);
            Global_String = strcspn((char *) buffer, "\"");

           /*  This loop is to strip out the colon in front of the 
               Variable Name */

            strncpy((char *) Variables[VarCount], " ", 319);
            for (j = 1; j < (int) strlen(STR.Strings[0]); j++) {
               Variables[VarCount][j-1] = STR.Strings[0][j]; 
            }
            VarCount += 1;
      
           /*  Now copy the Number */

            strcpy (Variables[VarCount], STR.Strings[1]); 
            VarCount += 1;
   
            if (STR.EqualSign) {
               OutSpaces[NumberOfLines] = Spaces;

               VariableMatch = 1;
               for (i = 0; i < NumberOfVariables; i++) {
                  if (strcmp (Variables[VarCount-2], OutVariables[i]) == 0) {
                     VariableMatch = 0;
                     break;
                  }
               }
               if (strcmp (Variables[VarCount-2], "history") == 0) {
                  VariableMatch = 0;
                  sprintf (OutLines[NumberOfLines], "HistoryCount = 1\n");
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
               }
               if (strcmp (Variables[VarCount-2], "title") == 0) {
                  VariableMatch = 0;
                  sprintf (OutLines[NumberOfLines], "TitleCount = 1\n");
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
               }
               if (VariableMatch) {
                  strcpy (OutVariables[NumberOfVariables], 
                          Variables[VarCount-2]);
                  NumberOfVariables += 1;
               }
               if (Global_String < Global_Length) {
                /*  The Global Variable is a string  */
                 if (VariableMatch) {
                    write_blanks (OutputFile, 7);
                    fprintf (OutputFile, "character*%d %s\n", 
                             strlen(Variables[VarCount-1]), 
                             Variables[VarCount-2]);
                 }
                 sprintf (OutLines[NumberOfLines], "%s = \"%s\"\n", 
                          Variables[VarCount-2], 
                          Variables[VarCount-1]);
                 OutSpaces[NumberOfLines] = Spaces;
                 NumberOfLines += 1;

               }
               else {
                /*  The Global Variable is an integer  */
                 if (VariableMatch) {
                    write_blanks (OutputFile, 7);
                    fprintf (OutputFile, "integer %s\n", Variables[VarCount-2]);
                 }
                 sprintf (OutLines[NumberOfLines], "%s = %s\n", 
                          Variables[VarCount-2], 
                          Variables[VarCount-1]);
                 OutSpaces[NumberOfLines] = Spaces;
                 NumberOfLines += 1;
               }
            }
            break;

         case 3:

           /*  Read in the definitions for the different variables */

            separate(buffer, &STR);

            for (i = 0; i < 50; i++) 
               strncpy((char *) SpecVar[SpecCount][i], " ", 319);

           /*  Copy it one character at a time and remove the 
               name in front of the variable */

            for (i = (int) strlen(SpecName[SpecCount-1])+1; 
                 i < (int) strlen(STR.Strings[0]); i++) {
               SpecVar[SpecCount][0][i-(strlen(SpecName[SpecCount-1])+1)] 
                      = STR.Strings[0][i];
            }

           /*  Copy the other strings */

            for (i = 1; i < STR.NumOfStrings; i++) {
               strcpy(SpecVar[SpecCount][i], STR.Strings[i]);
            }

            if (strcmp(SpecVar[SpecCount][0], "long_name") == 0 ||
                strcmp(SpecVar[SpecCount][0], "units") == 0     ||
                strcmp(SpecVar[SpecCount][0], "Geolocation_Pointer") == 0 ||
                strcmp(SpecVar[SpecCount][0], "Parameter_Type") == 0) {
             
             /*  Copy the string */

               if (strcmp(SpecVar[SpecCount][0], "units") == 0)
                  sprintf (OutLines[NumberOfLines],
                           "%s_str = \"%s\" \n", SpecVar[SpecCount][0],
                                           SpecVar[SpecCount][1]);
               else if (strcmp(SpecVar[SpecCount][0], "long_name") == 0) {

                  if (STR.LastChar == ',') {
                     sprintf (OutLines[NumberOfLines],
                              "%s = \n     $ \"%s\" //\n", SpecVar[SpecCount][0],
                                                         SpecVar[SpecCount][1]);
                     OutSpaces[NumberOfLines] = Spaces;

                     while (STR.LastChar == ',') {
                        getline(InputFile, buffer, 320);
                        separate(buffer, &STR);
                        NumberOfLines += 1;                  

                        if (STR.LastChar == ',') { 
                           sprintf (OutLines[NumberOfLines],
                                    "     $ \"%s\" //\n", STR.Strings[0]);         
                           OutSpaces[NumberOfLines] = 0;
                        }
                        else {
                           Spaces = 0;              
                           sprintf (OutLines[NumberOfLines],
                                    "     $ \"%s\"\n", STR.Strings[0]);                 
                        }
                     } 

                  }
                  else {
                     if (strlen(SpecVar[SpecCount][1]) > 90) {
                        strcpy (STR.Strings[0],SpecVar[SpecCount][1]);
                        trim_line(&STR, &NLH);
                        sprintf (SpecVar[SpecCount][0],
                              "%s = \n     $ \"%s\"", SpecVar[SpecCount][0],
                                                      NLH.Strings[0]);
                        for (i = 1; i < NLH.Lines; i++) {
                           sprintf (SpecVar[SpecCount][0],
                                   "%s //\n     $ \"%s\"", 
                                    SpecVar[SpecCount][0], NLH.Strings[i]);
                        }
                        sprintf (OutLines[NumberOfLines],"%s\n",
                                 SpecVar[SpecCount][0]);
                     }
                     else 
                        sprintf (OutLines[NumberOfLines],
                              "%s = \n     $ \"%s\" \n", SpecVar[SpecCount][0],
                                                         SpecVar[SpecCount][1]);
                  }
               }
               else {
                  if (strcmp(SpecVar[SpecCount][0], "Parameter_Type") == 0)
                     sprintf (OutLines[NumberOfLines],
                              "Param_Type = \"%s\" \n", SpecVar[SpecCount][1]);
                  else
                     sprintf (OutLines[NumberOfLines],
                              "%s = \"%s\" \n", SpecVar[SpecCount][0],
                                                SpecVar[SpecCount][1]);
               }
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
               Spaces = 13;
            }
            else if (strcmp(SpecVar[SpecCount][0], "valid_range") == 0 ||
                     strcmp(SpecVar[SpecCount][0], "_FillValue") == 0  ||
                     strcmp(SpecVar[SpecCount][0], "scale_factor") == 0 ||
                     strcmp(SpecVar[SpecCount][0], "add_offset") == 0) {
       
               if (SpecVar[SpecCount][1][0] == '\'') {
             
                /*  If the first character is a ' then copy the string */

                  strncpy((char *) hold, " ", 319);
                  for (i = 1; i < STR.NumOfStrings-1; i++) {
                     Length = strlen(SpecVar[SpecCount][i]);
                     Value = 0;

                    /*  Convert the first couple of strings, which is an 
                        octal number to a decimal number and write it out */

                     for (j = Length-3; j > 1; j--) {
                        hold[0] = SpecVar[SpecCount][i][j];
                        Value += atol(hold) 
                               * pow((double) 8, 
                                      fabs((double) (j - (Length - 3))));
                     }
                     if (strcmp(SpecVar[SpecCount][0], "_FillValue") == 0)
                        sprintf (OutLines[NumberOfLines],
                                 "Fill_Value = %d \n", Value);
                     else if (
                         strcmp(SpecVar[SpecCount][0], "scale_factor") == 0 ||
                         strcmp(SpecVar[SpecCount][0], "add_offset") == 0)
                        sprintf (OutLines[NumberOfLines],
                             "%s = %d \n", SpecVar[SpecCount][0], Value);
                     else
                        sprintf (OutLines[NumberOfLines],
                             "%s(%d) = %d \n", SpecVar[SpecCount][0], i, Value);
                     OutSpaces[NumberOfLines] = Spaces;
                     NumberOfLines += 1;
                     strncpy((char *) hold, " ", 319);
                  }

                 /*  Convert the last string, which is an 
                     octal number to a decimal number and write it out */

                  Length = strlen(SpecVar[SpecCount][STR.NumOfStrings-1]);
                  Value = 0;
                  for (j = Length-2; j > 1; j--) {
                     hold[0] = SpecVar[SpecCount][i][j];
                     Value += atol(hold)
                           * pow((double) 8, fabs((double) (j - (Length - 2))));
                  }
                  if (strcmp(SpecVar[SpecCount][0], "_FillValue") == 0)
                     sprintf (OutLines[NumberOfLines],
                              "Fill_Value = %d \n", Value);
                  else if (
                      strcmp(SpecVar[SpecCount][0], "scale_factor") == 0 ||
                      strcmp(SpecVar[SpecCount][0], "add_offset") == 0)
                     sprintf (OutLines[NumberOfLines],
                          "%s = %d \n", SpecVar[SpecCount][0], Value);
                  else
                     sprintf (OutLines[NumberOfLines],
                             "%s(%d) = %d \n", SpecVar[SpecCount][0], i, Value);
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
                  strncpy((char *) hold, " ", 319);
               }

              /*  
               *  If the data is a string, not a number then there is an
               *  error in the spec file (To Be Determined) was previously
               *  found 
               */

               else if (SpecVar[SpecCount][1][0] == 'T' ||
                        SpecVar[SpecCount][1][0] == 'p') {
                  for (i = 1; i < 3; i++) {
                     sprintf (OutLines[NumberOfLines],
                              "%s(%d) = %d \n", SpecVar[SpecCount][0], i, 0);
                     OutSpaces[NumberOfLines] = Spaces;
                     NumberOfLines += 1;
                  }
               }
               else {
                
                /*  Copy the string, which is a number and remove the last 
                    character and write it out*/
   
                 /*  The first one or two */

                  strncpy((char *) hold, " ", 319);
                  for (i = 1; i < STR.NumOfStrings-1; i++) {
                     Length = strlen(SpecVar[SpecCount][i]);
                     for (j = 0; j < Length-2; j++) {
                        hold[j] = SpecVar[SpecCount][i][j]; 
                     }
                     if (strcmp(SpecVar[SpecCount][0], "_FillValue") == 0)
                        sprintf (OutLines[NumberOfLines],
                                 "Fill_Value = %s \n", hold);
                     else if (
                         strcmp(SpecVar[SpecCount][0], "scale_factor") == 0 ||
                         strcmp(SpecVar[SpecCount][0], "add_offset") == 0)
                        sprintf (OutLines[NumberOfLines],
                             "%s = %s \n", SpecVar[SpecCount][0], hold);
                     else
                        sprintf (OutLines[NumberOfLines],
                              "%s(%d) = %s \n", SpecVar[SpecCount][0], i, hold);
                     OutSpaces[NumberOfLines] = Spaces;
                     NumberOfLines += 1;
                     strncpy((char *) hold, " ", 319);
                  }
   
                 /*  The last one */

                  strncpy((char *) hold, " ", 319);
                  Length = strlen(SpecVar[SpecCount][STR.NumOfStrings-1]);
                  for (j = 0; j < Length-1; j++) {
                     hold[j] = SpecVar[SpecCount][STR.NumOfStrings-1][j]; 
                  }
                  if (strcmp(SpecVar[SpecCount][0], "_FillValue") == 0)
                     sprintf (OutLines[NumberOfLines],
                              "Fill_Value = %s \n", hold);
                  else if (
                      strcmp(SpecVar[SpecCount][0], "scale_factor") == 0 ||
                      strcmp(SpecVar[SpecCount][0], "add_offset") == 0)
                     sprintf (OutLines[NumberOfLines],
                          "%s = %s \n", SpecVar[SpecCount][0], hold);
                  else
                     sprintf (OutLines[NumberOfLines],
                              "%s(%d) = %s \n", SpecVar[SpecCount][0], 
                               STR.NumOfStrings-1, hold);
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
                  strncpy((char *) hold, " ", 319);
               }
            }
            else if (strcmp(SpecVar[SpecCount][0], "Cell_Across_Swath_Sampling")
                     == 0 ||
                     strcmp(SpecVar[SpecCount][0], "Cell_Along_Swath_Sampling")
                     == 0) {

              /*  Store the numbers into the separate array spaces */

               strncpy((char *) hold, " ", 319);
               for (i = 1; i < STR.NumOfStrings-1; i++) {
                  Length = strlen(SpecVar[SpecCount][i]);

                 /*  Write out the first couple of numbers
                     (need to separate the , from the number */

                  for (j = 0; j < Length-1; j++) {
                     hold[j] = SpecVar[SpecCount][i][j];
                  }
                  sprintf (OutLines[NumberOfLines],
                           "%s(%d) = %s \n", SpecVar[SpecCount][0], i, hold);
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
                  strncpy((char *) hold, " ", 319);
               }

             /*  Write out the last number, no comma after this number */
   
               strncpy((char *) hold, " ", 319);
               Length = strlen(SpecVar[SpecCount][STR.NumOfStrings-1]);
               for (j = 0; j < Length; j++) {
                  hold[j] = SpecVar[SpecCount][STR.NumOfStrings-1][j];
               }
               sprintf (OutLines[NumberOfLines],
                        "%s(%d) = %s \n", SpecVar[SpecCount][0],
                         STR.NumOfStrings-1, hold);
               OutSpaces[NumberOfLines] = Spaces;
               NumberOfLines += 1;
               strncpy((char *) hold, " ", 319);
               NumberOfLines += 1;
            }
            else if (strcmp(SpecVar[SpecCount][0], "description") == 0 ||
                     strncmp((char *) SpecVar[SpecCount][0], "byte", 4) == 0 ||
                     DescCount > 0) {

               separate(buffer, &STR);

              /*  Write out the first description string */

               if (DescCount == 0) {
                  DescCount += 1;
                  sprintf (OutLines[NumberOfLines],
                       "description(%d) = \"%s\"\n", DescCount, STR.Strings[1]);
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
               }
               else if (STR.NumOfStrings > 0) {
 
                 /*  Call trim line to break the line into fortran usable 
                     lengths */

                  trim_line(&STR, &NLH);

                 /*  Write out each description string */

                  if (NLH.Lines > 1) {
                     DescCount += 1;
                     OutSpaces[NumberOfLines] = Spaces;
                     sprintf (OutLines[NumberOfLines],
                             "description(%d) = \n     $ \"%s\" // \n",
                              DescCount, NLH.Strings[0]);
                     for (i = 1; i < NLH.Lines-1; i++) {
                        sprintf (OutLines[NumberOfLines],
                          "%s     $ \"%s\"//\n", OutLines[NumberOfLines], 
                                                     NLH.Strings[i]);
                     }
                     sprintf (OutLines[NumberOfLines],
                             "%s     $ \"%s\"\n", OutLines[NumberOfLines], 
                                                  NLH.Strings[NLH.Lines-1]);
                     NumberOfLines += 1;
                  }
                  else {
                     DescCount += 1;
                     sprintf (OutLines[NumberOfLines],
                          "description(%d) = \n     $ \"%s\" \n", DescCount,
                           NLH.Strings[0]);
                     OutSpaces[NumberOfLines] = Spaces;
                     NumberOfLines += 1;
                  }
               }
            }
            break;
         
         case 4:

           /*  Read in the title */
   
            separate(buffer, &STR);

            if (STR.NumOfStrings > 0) {
               trim_line(&STR, &NLH);

              /*  Write out the title strings */

/*
               if (NLH.Lines > 1) {
                  TitleCount += 1;
                  OutSpaces[NumberOfLines] = Spaces;
                  sprintf (OutLines[NumberOfLines],
                          "title(%d) = \n     $ \"%s\" // \n",
                           TitleCount, NLH.Strings[0]);
                  for (i = 1; i < NLH.Lines-1; i++) {
                     sprintf (OutLines[NumberOfLines],
                       "%s     $ \"%s\" // \n", OutLines[NumberOfLines],
                                                  NLH.Strings[i]);
                  }
                  sprintf (OutLines[NumberOfLines],
                          "%s     $ \"%s\"\n", OutLines[NumberOfLines],
                                               NLH.Strings[NLH.Lines-1]);
                  NumberOfLines += 1;
               }
               else {
*/
                  TitleCount += 1;
                  sprintf (OutLines[NumberOfLines],
                       "title(%d) = \n     $ \"%s\" \n", TitleCount,
                        NLH.Strings[0]);
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
/*
               }
*/
            }
            break;
 
         case 5:

           /*  Read in the history */

            separate(buffer, &STR);

            if (STR.NumOfStrings > 0) {
               trim_line(&STR, &NLH);
  
              /*  Write the history strings out */

/*
               if (NLH.Lines > 1) {
                  HistoryCount += 1;
                  OutSpaces[NumberOfLines] = Spaces;
                  sprintf (OutLines[NumberOfLines],
                          "history(%d) = \n     $ \"%s\" // \n",
                           HistoryCount, NLH.Strings[0]);
                  for (i = 1; i < NLH.Lines-1; i++) {
                     sprintf (OutLines[NumberOfLines],
                       "%s     $ \"%s\" // \n", OutLines[NumberOfLines],
                                                  NLH.Strings[i]);
                  }
                  sprintf (OutLines[NumberOfLines],
                          "%s     $ \"%s\"\n", OutLines[NumberOfLines],
                                               NLH.Strings[NLH.Lines-1]);
                  NumberOfLines += 1;
               }
               else {
*/
                  HistoryCount = 1;
                  sprintf (OutLines[NumberOfLines],
                       "history(%d) = \n     $ \"%s\" \n", HistoryCount,
                        NLH.Strings[0]);
                  OutSpaces[NumberOfLines] = Spaces;
                  NumberOfLines += 1;
/*
               }
*/
            }
            break;

         case 6:

           /*  Read in the data list */

            separate(buffer, &STR);

           /*  Now copy the Number */

            strcpy (Variables[VarCount], STR.Strings[1]);
            DataCount = 0;
            DataVariables += 1;

            sprintf (Prog2Lines[Prog2Count], "   DataNames(%d) = \"%s\" ",
                     DataVariables, STR.Strings[0]);
            Prog2Count += 1;
            if (STR.EqualSign) {
               sprintf (Prog2Lines[Prog2Count], "   DataNames(%d) = %s",
                        DataVariables, STR.Strings[0]);
               for (i = 1; i < STR.NumOfStrings-1; i++) {
                  DataCount += 1;
                  Length = strlen(STR.Strings[i]);
                  strncpy((char *) hold, " ", 319);

                  for (j = 0; j < Length-1; j++) 
                     hold[j] = STR.Strings[i][j];

                  sprintf (Prog2Lines[Prog2Count], "   DataValues(%d,%d) = %s",
                           DataVariables, DataCount, hold);

                  Prog2Count += 1;
               }
               DataCount += 1;
               strcpy (hold, STR.Strings[STR.NumOfStrings-1]);
               sprintf (Prog2Lines[Prog2Count], "   DataValues(%d,%d) = %s",
                        DataVariables, DataCount, hold);

               Prog2Count += 1;
            }
            sprintf (Prog2Lines[Prog2Count], "   DataCount(%d) = %d",
                     DataVariables, DataCount);

            Prog2Count += 1;
            break;
         }
      }
        
     /*  Close off the fortran if statements */

      Spaces -= 3;
      OutSpaces[NumberOfLines] = Spaces;
      sprintf (OutLines[NumberOfLines], "End If\n");
      NumberOfLines += 1;

     /*  Write out the names of all the data type names */

      for (i = 0; i < SpecCount-1; i++) {
         sprintf (Prog2Lines[Prog2Count], "   FieldNames(%d) = \"%s\"", i+1, 
                  SpecName[i]);
         Prog2Count += 1;
      }
      sprintf (Prog2Lines[Prog2Count], "   FieldNames(%d) = \"%s\"\n", i+1, 
               SpecName[SpecCount-1]);
      Prog2Count += 1;
      sprintf (Prog2Lines[Prog2Count], "   FieldCount = %d\n", SpecCount);
      Prog2Count += 1;
      sprintf (Prog2Lines[Prog2Count], "End If\n");
      Prog2Count += 1;


      Spaces -= 3;
      OutSpaces[NumberOfLines] = Spaces;
      sprintf (OutLines[NumberOfLines], "End If\n\n");
      NumberOfLines += 1;

   }

  /*  Write the last lines for the fortran code */

   OutSpaces[NumberOfLines] = Spaces;
   sprintf (OutLines[NumberOfLines], "return\n");
   NumberOfLines += 1;

   OutSpaces[NumberOfLines] = Spaces;
   sprintf (OutLines[NumberOfLines], "end\n\n\n\n");
   NumberOfLines += 1;
    
  /*  Write the first fortran routine */

   fprintf (OutputFile, "\n");
   for (i = 0; i < NumberOfLines; i++) {
      write_blanks (OutputFile, OutSpaces[i]);
      fprintf (OutputFile, "%s", OutLines[i]);
   }

  /*  Write the second fortran routine */

   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "SUBROUTINE GET_SPEC_DATA(PROGRAMNAME, FIELDNAMES, DIMNAMES,\n");
   write_blanks (OutputFile, 5);
   fprintf (OutputFile, "$ ");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "DIMNUMBERS, DIMCOUNT, DATANAMES, DATACOUNT,\n");
   write_blanks (OutputFile, 5);
   fprintf (OutputFile, "$ ");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "DATAVALUES, NUMOFDATA, LIMIT, FIELDCOUNT)\n\n");

   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "implicit none\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "save\n\n");


   fprintf (OutputFile, "c---------------------------------------------------------------------\n");
   fprintf (OutputFile, "C!F77\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!DESCRIPTION:\n");
   fprintf (OutputFile, "c  Subroutine which extracts the correct file spec information\n");
   fprintf (OutputFile, "c  for the given product (MOD35, MOD07, or MOD06)\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!INPUT PARAMETERS:\n");
   fprintf (OutputFile, "c ProgramName   Ouput hdf file name\n");
   fprintf (OutputFile, "c Limit         Maximum size of Arrays\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!OUTPUT PARAMETERS:\n");
   fprintf (OutputFile, "c FieldNames    HDF Field Names for SDS\n");
   fprintf (OutputFile, "c FieldCount    Number of Dimension Values\n");
   fprintf (OutputFile, "c DimNames      Dimension Names for SDS\n");
   fprintf (OutputFile, "c DimNumbers    Dimension Values\n");
   fprintf (OutputFile, "c DimCount      Number of Dimension Values\n");
   fprintf (OutputFile, "c DataNames     Data Names for SDS\n");
   fprintf (OutputFile, "c DataNumbers   Data Values\n");
   fprintf (OutputFile, "c NumOfData     Number of Data Values\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!REVISION HISTORY:\n");
   fprintf (OutputFile, "c   This FORTRAN source file is created by spec2fort.c\n");
   fprintf (OutputFile, "c   spec2fort.c was designed by Walter.Wolf@ssec.wisc.edu\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!TEAM-UNIQUE HEADER:\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!REFERENCES AND CREDITS:\n");
   fprintf (OutputFile, "c\n");
   fprintf (OutputFile, "c!END\n");
   fprintf (OutputFile, "c---------------------------------------------------------------------\n\n");

   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "character*(*) ProgramName\n");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "character*(*) FieldNames(*), DimNames(*), DataNames(*)\n");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "integer Limit, FieldCount, NumOfData\n"); 
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "integer DimNumbers(*), DataCount(*), DimCount\n\n");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "real DataValues(Limit, Limit)\n");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "integer unlimited\n\n");

   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "NumOfData = 0\n");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "DimCount = 0\n");
   write_blanks (OutputFile, 7);
   fprintf (OutputFile, "FieldCount = 0\n\n");

   for (i = 0; i < Prog2Count; i++) {
      write_blanks (OutputFile, Spaces);
      fprintf (OutputFile, "%s\n", Prog2Lines[i]);
   }

   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "return\n");
   write_blanks (OutputFile, Spaces);
   fprintf (OutputFile, "end\n");

 /*  Close the Files */

   fclose (InputFile);
   fclose (OutputFile);
   return(1);
}
         
   
