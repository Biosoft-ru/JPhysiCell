package ru.biosoft.physicell.core;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Map;

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

public class Rules
{
    static Map<CellDefinition, HypothesisRuleset> hypothesis_rulesets;

    static void add_hypothesis_ruleset(CellDefinition pCD)
    {

        //        auto search = hypothesis_rulesets.find( pCD );
        //        if( search == hypothesis_rulesets.end() )
        //        {
        HypothesisRuleset HRS = new HypothesisRuleset();
        HRS.sync_to_CellDefinition( pCD );
        //            hypothesis_rulesets[pCD] = HRS; 
        //        }
        hypothesis_rulesets.put( pCD, HRS );
    }

    static void intialize_hypothesis_rulesets()
    {
        hypothesis_rulesets.clear(); // empty(); 

        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            //            CellDefinition pCD = CellDefinitions_by_index[n];
            add_hypothesis_ruleset( cd );
        }
    }

    HypothesisRuleset access_ruleset(CellDefinition pCD)
    {
        return hypothesis_rulesets.get( pCD );
    }

    static HypothesisRuleset find_ruleset(CellDefinition pCD)
    {
        return ( hypothesis_rulesets.get( pCD ) );
    }

    static String display_hypothesis_rulesets()
    {
        StringBuilder sb = new StringBuilder();
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            sb.append( hypothesis_rulesets.get( cd ) );//CellDefinitions_by_index[n]] );
        }
        return sb.toString();
    }

    String detailed_display_hypothesis_rulesets()
    {
        StringBuilder sb = new StringBuilder();
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            sb.append( hypothesis_rulesets.get( cd ).detailed_display() );
        }
        return sb.toString();
    }

    void add_rule(String cell_type, String signal, String behavior, String response) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            System.out.println( "Warning: Attempted to add rule for " + cell_type + ", but no cell definition found for this type." );
        }

        HypothesisRuleset pHRS = find_ruleset( pCD );
        if( pHRS == null )
        {
            System.out.println( "Warning: Attempted to add rule for " + cell_type + ", but no hypothesis ruleset found for this type." );
        }

        if( pHRS.find_behavior( behavior ) != null )
        {
            if( ( pHRS ).get( behavior ).behavior != behavior )
            {
                ( pHRS ).get( behavior ).behavior = behavior;
                System.out.println( "wha?" );
            }
        }
        pHRS.add_behavior( behavior );
        ( pHRS ).get( behavior ).add_signal( signal, response );
    }

    void add_rule(String cell_type, String signal, String behavior, String response, boolean use_for_dead) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            throw new Exception( "Warning: Attempted to add rule for " + cell_type + ", but no cell definition found for this type." );
        }

        HypothesisRuleset pHRS = find_ruleset( pCD );
        if( pHRS == null )
        {
            throw new Exception( "Warning: Attempted to add rule for " + cell_type + ", but no hypothesis ruleset found for this type." );
        }

        if( pHRS.find_behavior( behavior ) != null )
        {
            if( ( pHRS ).get( behavior ).behavior != behavior )
            {
                ( pHRS ).get( behavior ).behavior = behavior;
                System.out.println( "wha?" );
            }
        }

        pHRS.add_behavior( behavior );

        pHRS.get( behavior ).add_signal( signal, response );

        // set dead flag

        int n = ( pHRS ).get( behavior ).find_signal( signal );
        pHRS.get( behavior ).applies_to_dead_cells.set( n, use_for_dead );
    }

    void set_hypothesis_parameters(String cell_type, String signal, String behavior, double half_max, double hill_power) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + cell_type
                    + ", but no cell definition found for this type." );
        }

        if( find_ruleset( pCD ) == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + cell_type
                    + ", but no behavior ruleset for this cell type." );
        }

        if( hypothesis_rulesets.get( pCD ).find_behavior( behavior ) == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + cell_type
                    + ", but the cell type has no rules for this behavior." );

            // System.out.println( "Error. No " + behavior + " rules for " + cell_type + ". Ignoring."+ std::endl; return;
        }

        if( hypothesis_rulesets.get( pCD ).get( behavior ).find_signal( signal ) < 0 )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + cell_type
                    + ", but the cell type's behavior does not vary with this signal." );
        }

        hypothesis_rulesets.get( pCD ).get( behavior ).set_half_max( signal, half_max );
        hypothesis_rulesets.get( pCD ).get( behavior ).set_hill_power( signal, hill_power );
    }

    void set_behavior_parameters(String cell_type, String behavior, double min_value, double max_value) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + cell_type
                    + ", but the cell definition is not found." );
        }

        if( find_ruleset( pCD ) == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + cell_type
                    + ", but there is no hypothesis ruleset for this cell type." );
        }

        if( hypothesis_rulesets.get( pCD ).find_behavior( behavior ) == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + cell_type
                    + ", but there is no rules for this behavior for this cell type." );
        }

        hypothesis_rulesets.get( pCD ).get( behavior ).min_value = min_value;
        hypothesis_rulesets.get( pCD ).get( behavior ).max_value = max_value;
    }

    void set_behavior_parameters(String cell_type, String behavior, double min_value, double base_value, double max_value) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + cell_type
                    + ", but the cell definition is not found." );
        }

        if( find_ruleset( pCD ) == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + cell_type
                    + ", but there is no hypothesis ruleset for this cell type." );
        }


        if( hypothesis_rulesets.get( pCD ).find_behavior( behavior ) == null )
        {
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + cell_type
                    + ", but there is no rules for this behavior for this cell type." );
        }

        hypothesis_rulesets.get( pCD ).get( behavior ).min_value = min_value;
        hypothesis_rulesets.get( pCD ).get( behavior ).max_value = max_value;
        hypothesis_rulesets.get( pCD ).get( behavior ).base_value = base_value;
    }


    void set_behavior_base_value(String cell_type, String behavior, double base_value) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but the cell definition is not found." );
        }

        if( find_ruleset( pCD ) == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but no hypothesis ruleset found for this cell type." );
        }

        if( hypothesis_rulesets.get( pCD ).find_behavior( behavior ) == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but no rule for this behavior found for this cell type." );
        }

        hypothesis_rulesets.get( pCD ).get( behavior ).base_value = base_value;
    }

    void set_behavior_min_value(String cell_type, String behavior, double min_value) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but the cell definition is not found." );
        }

        if( find_ruleset( pCD ) == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but no hypothesis ruleset found for this cell type." );
        }

        if( hypothesis_rulesets.get( pCD ).find_behavior( behavior ) == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but no rule for this behavior found for this cell type." );
        }

        hypothesis_rulesets.get( pCD ).get( behavior ).min_value = min_value;
    }

    void set_behavior_max_value(String cell_type, String behavior, double max_value) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        if( pCD == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but the cell definition is not found." );
        }

        if( find_ruleset( pCD ) == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but no hypothesis ruleset found for this cell type." );
        }

        if( hypothesis_rulesets.get( pCD ).find_behavior( behavior ) == null )
        {
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + cell_type
                    + ", but no rule for this behavior found for this cell type." );
        }

        hypothesis_rulesets.get( pCD ).get( behavior ).max_value = max_value;
    }

    void apply_ruleset(Cell pCell) throws Exception
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( pCell.type_name );
        hypothesis_rulesets.get( pCD ).apply( pCell );
    }

    void rule_phenotype_function(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        apply_ruleset( pCell );

        // by default, rules only apply to live cells
        // s

        /*
        // safety checks for dead cells 
        if( get_single_signal(pCell,"dead") > 0.11 )
        {
            // can't die twice
            set_single_behavior(pCell,"apoptosis",0.0);
            set_single_behavior(pCell,"necrosis",0.0);
        
            // can't cycle 
            set_single_behavior(pCell,"cycle entry",0.0);
        
            // can't crawl 
            set_single_behavior(pCell,"migration speed",0.0);
        }
        */
    }

    /* add these to core */
    /*
    double[] linear_response_to_Hill_parameters( double s0, double s1 )
    {
        static double tol = 0.1; 
        static double param1 = (1-tol)/tol; 
        static double param2 = log(param1); 
    
        // half max, then hill power 
        double hm = 0.5* (s0+s1); 
    
        // hp so that H(s1) ~ (1-tol)
        double hp = round( param2 / log(s1/hm) ); 
    
        double[] output = { hm , hp }; 
    
        return output; 
    }
    
    double[] Hill_response_to_linear_parameters( double half_max , double Hill_power )
    {
        static double tol = 0.1; 
        static double param1 = (1-tol)/tol; 
        double param2 = pow( param1 , 1.0/ Hill_power ); 
    
        // s1 such that H(s1) ~ (1-tol)
        double s1 = half_max * param2; 
    
        // s0 for symmetry
        double s0 = 2*half_max -s1; 
        if( s0 < 0 )
        { s0 = 0.0; }
    
        double[] output = {s0,s1}; 
    
        return output; 
    }
    */

    //    void split_csv( String input , String[] output , char delim )
    //    {
    //        input.split( delim );
    ////        output.resize(0); 
    ////
    ////        std::istringstream is(input);
    ////        String part;
    ////        while( getline(is, part, delim ) )
    ////        { output.push_back(part); }
    //    }

    String csv_strings_to_English_v1(String[] strings, boolean include_cell_header)
    {
        String output = "";

        if( include_cell_header )
        {
            output += "In ";
            output += strings[0];
            output += " cells:\n\t"; // In {cell type X} cells: 
        }

        // malignant epithelial,oxygen,decreases,necrosis,2.80E-03,0,decreases,3.75,8,0

        output += strings[1]; // {signal}
        output += " ";

        output += strings[2]; // {increases/decreases}
        output += " ";

        output += strings[3]; // {behavior}


        output += " from "; // {base}
        output += strings[4];

        output += " towards ";
        output += strings[5];

        output += " with a Hill response, with half-max ";
        output += strings[6];

        output += " and Hill power ";
        output += strings[7];

        output += ".";
        boolean use_when_dead = false;
        char start_char = Character.toUpperCase( strings[8].charAt( 0 ) );
        //        char start_char = toupper( strings[8][0] );
        if( start_char == 'T' || start_char == '1' )
        {
            output += " Rule applies to dead cells.";
        }

        return output;
    }

    String csv_strings_to_English_v2(String[] strings, boolean include_cell_header)
    {
        String output = "";

        if( include_cell_header )
        {
            output += "In ";
            output += strings[0];
            output += " cells:\n\t"; // In {cell type X} cells: 
        }

        output += strings[1]; // {signal}
        output += " ";

        output += strings[2]; // {increases/decreases}
        output += " ";

        output += strings[3]; // {behavior}


        //  output += " from "; // {base}
        //  output += strings[4]; 

        output += " towards ";
        output += strings[4];

        output += " with a Hill response, with half-max ";
        output += strings[5];

        output += " and Hill power ";
        output += strings[6];

        output += ".";
        boolean use_when_dead = false;
        char start_char = Character.toUpperCase( strings[7].charAt( 0 ) );//toupper( strings[7][0] );
        if( start_char == 'T' || start_char == '1' )
        {
            output += " Rule applies to dead cells.";
        }

        return output;
    }


    String csv_strings_to_English(String[] strings, boolean include_cell_header)
    {
        String output = "";

        if( include_cell_header )
        {
            output += "In ";
            output += strings[0];
            output += " cells:\n\t"; // In {cell type X} cells: 
        }

        output += strings[5]; // {signal}
        output += " ";

        output += strings[6]; // {increases/decreases}
        output += " ";

        output += strings[1]; // {behavior}


        output += " from "; // {base}
        output += strings[3];

        output += " towards ";
        if( strings[6].startsWith( "i" ) || strings[6].startsWith( "I" ) )
        {
            output += strings[4];
        }
        else
        {
            output += strings[2];
        }

        output += " with a Hill response, with half-max ";
        output += strings[7];

        output += " and Hill power ";
        output += strings[8];

        output += ".";
        boolean use_when_dead = false;
        char start_char = strings[9].charAt( 0 );// toupper( strings[9][0] );
        if( start_char == 'T' || start_char == '1' )
        {
            output += " Rule applies to dead cells.";
        }

        return output;
    }


    String csv_strings_to_English_HTML(String[] strings, boolean include_cell_header)
    {
        String output = "<p>";

        if( include_cell_header )
        {
            output += "In ";
            output += strings[0];
            output += " cells:<br>\n"; // In {cell type X} cells: 
        }

        output += "&nbsp;";
        output += strings[5]; // {signal}
        output += " ";

        output += strings[6]; // {increases/decreases}
        output += " ";

        output += strings[1]; // {behavior}

        output += " from "; // {base}
        output += strings[3];

        output += " towards ";
        if( strings[6].startsWith( "i" ) || strings[6].startsWith( "I" ) )
        {
            output += strings[4];
        }
        else
        {
            output += strings[2];
        }

        output += " with a Hill response, with half-max ";
        output += strings[7];

        output += " and Hill power ";
        output += strings[8];

        output += ".";
        boolean use_when_dead = false;
        char start_char = Character.toUpperCase( strings[9].charAt( 0 ) );//toupper( strings[9][0] );
        if( start_char == 'T' || start_char == '1' )
        {
            output += " Rule applies to dead cells.";
        }
        output += "\n</p>\n";

        return output;
    }


    /*
    0              1          2        3           4         5       6          7         8           9 
     Cell type, behavior, min value, base value, max value, signal, direction, half-max, Hill power, dead
    */

    void parse_csv_rule_v0(String[] input) throws Exception
    {
        // if it's wrong length, skip 
        boolean skip = false;
        if( input.length != 10 )
        {
            skip = true;
        }
        // if any empty strings, skip
        for( int n = 0; n < input.length; n++ )
        {
            if( input[n].isEmpty() == true )
            {
                skip = true;
            }
        }
        if( skip == true )
        {
            System.out.println( "Warning: Misformed rule (likely from an empty rules file). Skipping." );
            return;
        }

        String temp = csv_strings_to_English( input, false );

        // string portions of the rule
        String cell_type = input[0];
        String behavior = input[1];
        String signal = input[5];
        String response = input[6];

        // numeric portions of the rule 
        double min_value = Double.parseDouble( input[2] );//std::atof( input[2].c_str() );
        double base_value = Double.parseDouble( input[3] );//= std::atof( input[3].c_str() );
        double max_value = Double.parseDouble( input[4] );//= std::atof( input[4].c_str() ); 

        double half_max = Double.parseDouble( input[7] );//std::atof( input[7].c_str() );
        double hill_power = Double.parseDouble( input[8] );//std::atof( input[8].c_str() );
        boolean use_for_dead = Boolean.parseBoolean( input[9] );//(bool) std::atof( input[9].c_str() ); 

        System.out.println( "Adding rule for " + cell_type + " cells:\n\t" );
        System.out.println( temp );

        // add_rule(cell_type,signal,behavior,response);  
        add_rule( cell_type, signal, behavior, response, use_for_dead );

        set_hypothesis_parameters( cell_type, signal, behavior, half_max, hill_power );

        // compare to base behavior value in cell def for discrepancies 
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        double ref_base_value = SignalBehavior.get_single_base_behavior( pCD, behavior );
        if( Math.abs( ref_base_value - base_value ) > 1e-15 )
        {
            throw new Exception( "Error: Base value for " + behavior + " in cell type " + cell_type//  + std::endl 
                    + "       has base value " + base_value + " in the rule, " // + std::endl 
                    + "but base value " + ref_base_value + " in the cell definition." //+ std::endl
                    + "       Fix this discrepancy to continue the model." );
        }

        // get_single_base_behavior

        set_behavior_parameters( cell_type, behavior, min_value, base_value, max_value );

        /*
        // set "applies to dead"
        
        int n = access_ruleset(pCD)[behavior].find_signal(signal); 
        access_ruleset(pCD)[behavior].applies_to_dead_cells[n] = use_for_dead; 
        */
    }

    void parse_csv_rule_v0(String input) throws Exception
    {
        String[] tokenized_string = input.split( "," );
        //        split_csv( input, tokenized_string, ',' );

        // Make sure it was truly comma-separated. If not, try tab.
        if( tokenized_string.length != 10 )
        {
            tokenized_string = input.split( "\t" );
            //            split_csv( input, tokenized_string, '\t' );
        }

        parse_csv_rule_v0( tokenized_string );
    }

    void parse_csv_rules_v0(String filename) throws Exception
    {
        //        std::fstream fs( filename, std::ios::in );
        File f = new File( filename );
        if( !f.exists() )
        {
            System.out.println( "Warning: Rules file " + filename + " failed to open." );
            return;
        }

        System.out.println( "Processing rules in file " + filename + " ... " );

        try (BufferedReader br = new BufferedReader( new FileReader( f ) ))
        {
            String line = br.readLine();
            while( line != null )
            {
                if( line.length() > 0 )
                {
                    parse_csv_rule_v0( line );
                }
            }
            line = br.readLine();

        }

        //        }while(fs.eof()==false){
        //            String line;
        //        std::getline(fs,line,'\n');
        //        if(line.size()>0){parse_csv_rule_v0(line);}}

        //        fs.close();

        System.out.println( "Done!" );
    }

    /* v1 work */

    /*
     v0:::
    0              1          2        3           4         5       6          7         8           9 
     Cell type, behavior, min value, base value, max value, signal, direction, half-max, Hill power, dead
    
    
     v1:::
     0          1       2          3         4           5                   6         7           8
     Cell type, signal, direction, behavior, base value, max response value, half-max, Hill power, applies to dead? 
    */

    void parse_csv_rule_v1(String[] input) throws Exception
    {
        // if it's wrong length, skip 
        boolean skip = false;
        if( input.length != 9 )
        {
            skip = true;
        }
        // if any empty strings, skip
        for( int n = 0; n < input.length; n++ )
        {
            if( input[n].isEmpty() == true )
            {
                skip = true;
            }
        }
        if( skip == true )
        {
            System.out.println( "Warning: Misformed rule (likely from an empty rules file). Skipping." );

            for( int n = 0; n < input.length; n++ )
            {
                System.out.println( n + " : " + input[n] );
            }

            return;
        }

        String temp = csv_strings_to_English_v1( input, false ); // need a v1 version of this

        // string portions of the rule
        String cell_type = input[0];
        String signal = input[1];
        String response = input[2];
        String behavior = input[3];

        // numeric portions of the rule 
        // double min_value  = std::atof( input[2].c_str() );

        double base_value = Double.parseDouble( input[4] );//std::atof( input[4].c_str() );
        double max_response = Double.parseDouble( input[5] );//std::atof( input[5].c_str() ); 

        // hmm from here 
        // double max_value  = std::atof( input[4].c_str() ); 

        double half_max = Double.parseDouble( input[6] );//  std::atof( input[6].c_str() );
        double hill_power = Double.parseDouble( input[7] );//std::atof( input[7].c_str() );
        boolean use_for_dead = Boolean.parseBoolean( input[8] );//(bool) std::atof( input[8].c_str() ); 

        System.out.println( "Adding rule for " + cell_type + " cells:\n\t" );
        System.out.println( temp );

        // add_rule(cell_type,signal,behavior,response);  
        add_rule( cell_type, signal, behavior, response, use_for_dead );
        set_hypothesis_parameters( cell_type, signal, behavior, half_max, hill_power );

        // compare to base behavior value in cell def for discrepancies 
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        double ref_base_value = SignalBehavior.get_single_base_behavior( pCD, behavior );
        if( Math.abs( ref_base_value - base_value ) > 1e-15 )
        {
            throw new Exception( "Error: Base value for " + behavior + " in cell type " + cell_type + "\n" + "       has base value "
                    + base_value + " in the rule, " + "\n"// + std::endl 
                    + "but base value " + ref_base_value + " in the cell definition." + "\n"
                    + "       Fix this discrepancy to continue the model." );
        }

        set_behavior_base_value( cell_type, behavior, base_value );

        if( response == "increases" )
        {
            set_behavior_max_value( cell_type, behavior, max_response );
        }
        else
        {
            set_behavior_min_value( cell_type, behavior, max_response );
        }
    }

    void parse_csv_rule_v1(String input) throws Exception
    {
        String[] tokenized_string = input.split( "," );
        //        split_csv( input, tokenized_string, ',' );

        // Make sure it was truly comma-separated. 
        // If not, try tab.
        if( tokenized_string.length != 9 )
        {
            tokenized_string = input.split( ",\t" );
            //            split_csv( input, tokenized_string, '\t' );
        }

        // check for comment 
        if( tokenized_string[0].charAt( 0 ) == '/' && tokenized_string[0].charAt( 1 ) == '/' )
        {
            System.out.println( "Skipping commented rule (" + input + ")" );
            return;
        }

        parse_csv_rule_v1( tokenized_string );
    }

    void parse_csv_rules_v1(String filename) throws Exception
    {
        //        std::fstream fs( filename, std::ios::in );
        File f = new File( filename );
        if( !f.exists() )
        {
            System.out.println( "Warning: Rules file " + filename + " failed to open." );
            return;
        }

        System.out.println( "Processing rules in file " + filename + " ... " );

        try (BufferedReader br = new BufferedReader( new FileReader( f ) ))
        {
            String line = br.readLine();
            while( line != null )
            {
                if( line.length() > 0 )
                    parse_csv_rule_v1( line );
                line = br.readLine();
            }
        }
        //        while( fs.eof() == false )
        //        {
        //            String line;   
        //            std::getline( fs , line, '\n'); 
        //            if( line.size() > 0 )
        //            { parse_csv_rule_v1(line); }
        //        }
        //
        //        fs.close(); 

        System.out.println( "Done!" );
    }


    /* end of v1 work */

    /* v2 work */

    /*
     v0:::
    0              1          2        3           4         5       6          7         8           9 
     Cell type, behavior, min value, base value, max value, signal, direction, half-max, Hill power, dead
    
    
     v1:::
     0          1       2          3         4           5                   6         7           8
     Cell type, signal, direction, behavior, base value, max response value, half-max, Hill power, applies to dead? 
    
     v2:::
     0          1       2          3         4                   5         6           7           
     Cell type, signal, direction, behavior, max response value, half-max, Hill power, applies to dead?  
    */

    void parse_csv_rule_v2(String[] input) throws Exception
    {
        // if it's wrong length, skip 
        boolean skip = false;
        if( input.length != 8 )
        {
            skip = true;
        }
        // if any empty strings, skip
        for( int n = 0; n < input.length; n++ )
        {
            if( input[n].isEmpty() == true )
            {
                skip = true;
            }
        }
        if( skip == true )
        {
            System.out.println( "Warning: Misformed rule (likely from an empty rules file). Skipping." );

            for( int n = 0; n < input.length; n++ )
            {
                System.out.println( n + " : " + input[n] );
            }
            return;
        }

        String temp = csv_strings_to_English_v2( input, false ); // need a v1 version of this

        // string portions of the rule
        String cell_type = input[0];
        String signal = input[1];
        String response = input[2];
        String behavior = input[3];

        // numeric portions of the rule 
        // double min_value  = std::atof( input[2].c_str() );

        // double base_value = std::atof( input[4].c_str() );
        double max_response = Double.parseDouble( input[4] );// .std::atof(input[4].c_str());

        // hmm from here 
        // double max_value  = std::atof( input[4].c_str() ); 

        double half_max = Double.parseDouble( input[5] );
        double hill_power = Double.parseDouble( input[6] );
        boolean use_for_dead = Boolean.parseBoolean( input[7] );//::atof(input[7].c_str());

        System.out.println( "Adding rule for " + cell_type + " cells:\n\t" );
        System.out.println( temp );

        // add_rule(cell_type,signal,behavior,response);  
        add_rule( cell_type, signal, behavior, response, use_for_dead );
        set_hypothesis_parameters( cell_type, signal, behavior, half_max, hill_power );

        // compare to base behavior value in cell def for discrepancies 
        CellDefinition pCD = CellDefinition.getCellDefinition( cell_type );
        double ref_base_value = SignalBehavior.get_single_base_behavior( pCD, behavior );

        set_behavior_base_value( cell_type, behavior, ref_base_value );
        if( response == "increases" )
        {
            set_behavior_max_value( cell_type, behavior, max_response );
        }
        else
        {
            set_behavior_min_value( cell_type, behavior, max_response );
        }
    }

    void parse_csv_rule_v2(String input) throws Exception
    {
        String[] tokenized_string = input.split( "," );
        //        split_csv(input,tokenized_string,',');

        // Make sure it was truly comma-separated. 
        // If not, try tab.
        if( tokenized_string.length != 8 )
        {
            tokenized_string = input.split( "\t" );
        }

        // check for comment 
        if( tokenized_string[0].charAt( 0 ) == '/' && tokenized_string[0].charAt( 1 ) == '/' )
        {
            System.out.println( "Skipping commented rule (" + input + ")" );
            return;
        }

        parse_csv_rule_v2( tokenized_string );
    }

    void parse_csv_rules_v2(String filename) throws Exception
    {
        File f = new File( filename );
        if( !f.exists() )
        {
            System.out.println( "Warning: Rules file " + filename + " failed to open." );
            return;
        }

        System.out.println( "Processing rules in file " + filename + " ... " );

        try (BufferedReader br = new BufferedReader( new FileReader( f ) ))
        {
            String line = br.readLine();
            while( line != null )
            {
                if( line.length() > 0 )
                {
                    parse_csv_rule_v2( line );
                }
                line = br.readLine();
            }
        }
        //        while( fs.eof() == false )
        //        {
        //            String line;   
        //            std::getline( fs , line, '\n'); 
        //            if( line.size() > 0 )
        //            { parse_csv_rule_v2(line); }
        //        }
        //
        //        fs.close(); 

        System.out.println( "Done!" );
    }


    /* end of v2 work */

    // needs fixing
    //    void parse_rules_from_pugixml()
    //    {
    //        pugi::xml_node node=physicell_config_root.child("cell_rules");if(!node){System.out.println("Error: Could not find <cell_rules> section of XML config file.\n"+"       Cannot parse cell rules, so disabling.");
    //
    //        PhysiCell_settings.rules_enabled=false;return;}
    //
    //        // find the set of rulesets 
    //        node=node.child("rulesets");if(!node){System.out.println("Error: Could not find <rulesets> in the <cell_rules> section of XML config file.\n"+"       Cannot parse cell rules, so disabling.");
    //
    //        PhysiCell_settings.rules_enabled=false;return;}
    //        // find the first ruleset 
    //        node=node.child("ruleset");if(!node){System.out.println("Error: Could not find any <ruleset> in the <rulesets> section of XML config file.\n"+"       Cannot parse cell rules, so disabling.");
    //
    //        PhysiCell_settings.rules_enabled=false;return;}
    //
    //        while(node){System.out.println(node.name());if(node.attribute("enabled").as_bool()==true){String folder=xml_get_string_value(node,"folder");String filename=xml_get_string_value(node,"filename");String input_filename=folder+"/"+filename;
    //
    //        System.out.println("\tProcessing ruleset in "+input_filename+" ... ");String format=node.attribute("format").as_string();String protocol=node.attribute("protocol").as_string();double version=node.attribute("version").as_double();
    //
    //        boolean done=false;
    //
    //        if(format=="CSV"||format=="csv"){if(version<1.0){System.out.println("\tFormat: CSV (prototype version)");
    //
    //        parse_csv_rules_v0(input_filename); // parse all rules in a CSV file 
    //
    //        PhysiCell_settings.rules_enabled=true;
    //
    //        done=true;}
    //
    //        if(version>=1.0-1e-10&&version<2.0-1e-10&&protocol=="CBHG"&&done==false){System.out.println("\tFormat: CSV (version "+version+")");
    //
    //        parse_csv_rules_v1(input_filename); // parse all rules in a CSV file 
    //
    //        PhysiCell_settings.rules_enabled=true;
    //
    //        done=true;}
    //
    //        if(version>=2.0-1e-10&&protocol=="CBHG"&&done==false){System.out.println("\tFormat: CSV (version "+version+")");
    //
    //        parse_csv_rules_v2(input_filename); // parse all rules in a CSV file 
    //
    //        PhysiCell_settings.rules_enabled=true;
    //
    //        done=true;}
    //
    //        }
    //
    //
    //        if(done==false){System.out.println("\tWarning: Ruleset had unknown format ("+format+"). Skipping!");}
    //
    //        }else{System.out.println("\tRuleset disabled ... ");}node=node.next_sibling("ruleset");}return;
    //
    //        exit(0);
    //
    //        // enabled? 
    //        if(node.attribute("enabled").as_bool()==false){return;}
    //
    //        // get filename 
    //
    //        String folder=xml_get_string_value(node,"folder");String filename=xml_get_string_value(node,"filename");String input_filename=folder+"/"+filename;
    //
    //        String filetype=node.attribute("type").value();
    //
    //        // what kind? 
    //        if(filetype=="csv"||filetype=="CSV"){System.out.println("Loading rules from CSV file "+input_filename+" ... ");
    //        // load_cells_csv( input_filename );
    //        parse_csv_rules_v0(input_filename);return;}
    //
    //        return;
    //    }

    void parse_rules_from_parameters_v0(Model model) throws Exception
    {
        boolean enabled = model.getParameterBoolean( "rules_enabled" );

        // enabled? 
        if( enabled == false )
        {
            return;
        }

        // get filename 

        String folder = model.getParameter( "rules_folder" );
        String filename = model.getParameter( "rules_filename" );
        String input_filename = folder + "/" + filename;

        String filetype = model.getParameter( "rules_type" );
        ;

        // what kind? 
        if( filetype == "csv" || filetype == "CSV" )
        {
            System.out.println( "Loading rules from CSV file " + input_filename + " ... " );
            // load_cells_csv( input_filename );
            parse_csv_rules_v0( input_filename );
            return;
        }
    }

    void spaces_to_underscore(String str)
    {
        str.replace( ' ', '_' );
        //        for( int n=0 ; n < str.size(); n++ )
        //        { if( str[n] == ' ' ){ str[n] = '_'; } }
    }

    /*
    submit this as bugfix in PhysiCell (PhysiCell_settings.cpp)
    
    template <class T>
    int Parameters<T>::find_index( String search_name )
    {
        auto search = name_to_index_map.find( search_name );
        if( search != name_to_index_map.end() )
        { return search.second; }
        return -1; 
        // return name_to_index_map[ search_name ]; 
    }
    
    */

    public static String stream_annotated_English_rules()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Cell Hypothesis Rules\n\n" );
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            HypothesisRuleset pHRS = find_ruleset( cd );
            sb.append( "In " + pHRS.cell_type + " cells:\n" );

            for( int k = 0; k < pHRS.rules.size(); k++ )
            {
                sb.append( pHRS.rules.get( k ).English_display() );
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }

    public static String stream_annotated_English_rules_HTML()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "<html>\n<body><h1>Cell Hypothesis Rules</h1>\n" );
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            sb.append( "<p>" );
            //        CellDefinition pCD=CellDefinitions_by_index[n];
            HypothesisRuleset pHRS = find_ruleset( cd );
            sb.append( "In " + pHRS.cell_type + " cells:\n" );
            sb.append( "<ul>\n" );
            for( int k = 0; k < pHRS.rules.size(); k++ )
            {
                sb.append( pHRS.rules.get( k ).English_display_HTML() );
            }
            sb.append( "</ul>\n</p>\n" );
        }
        sb.append( "</body>\n</html>\n" );
        return sb.toString();
    }

    static void save_annotated_English_rules(String filename)
    {
        //        String filename = PhysiCellSettings.folder + "/rules.txt";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_English_rules();
        //        of.close(); 
    }

    public static String stream_annotated_detailed_English_rules()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Cell Hypothesis Rules (detailed)\n\n" );
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            HypothesisRuleset pHRS = find_ruleset( cd );
            sb.append( "In " + pHRS.cell_type + " cells:\n" );
            for( int k = 0; k < pHRS.rules.size(); k++ )
            {
                pHRS.rules.get( k ).English_detailed_display();
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }

    static String stream_annotated_detailed_English_rules_HTML()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "<html>\n<body><h1>Cell Hypothesis Rules (detailed)</h1>\n" );// + std::endl; 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            sb.append( "<p>" );
            //            CellDefinition pCD = CellDefinitions_by_index[n]; 
            HypothesisRuleset pHRS = find_ruleset( cd );
            sb.append( "In " + pHRS.cell_type + " cells:\n" );// + std::endl; 
            sb.append( "<ul>" );// + std::endl; 
            for( int k = 0; k < pHRS.rules.size(); k++ )
            {
                sb.append( pHRS.rules.get( k ).English_detailed_display_HTML() );
            }
            sb.append( "</ul>\n</p>\n" );// + std::endl; 
        }
        sb.append( "</body>\n</html>\n" );// + std::endl; 
        return sb.toString();
    }

    static void save_annotated_detailed_English_rules(String filename) throws Exception
    {
        //        String filename = PhysiCell_settings.folder + "/detailed_rules.txt";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_detailed_English_rules();
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            bw.write( text );
        }
        //        of.close(); 
    }

    static void save_annotated_detailed_English_rules_HTML(String filename) throws Exception
    {
        //        String filename = PhysiCell_settings.folder + "/detailed_rules.html";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_detailed_English_rules_HTML();
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            bw.write( text );
        }
        //        of.close(); 
    }

    static void save_annotated_English_rules_HTML(String filename) throws Exception
    {
        //        String filename = PhysiCell_settings.folder + "/rules.html";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_English_rules_HTML();
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            bw.write( text );
        }
        //        of.close(); 
    }

    // v0 version 

    static void export_rules_csv_v0(String filename) throws Exception
    {
        File f = new File( filename );
        if( !f.exists() )
        {
            System.out.println( "Warning: Rules export file " + filename + " failed to open." );
            return;
        }

        System.out.println( "Exporting rules to file " + filename + " (v0 format) ... " );

        for( CellDefinition pCD : CellDefinition.getCellDefinitions() )// n=0; n < CellDefinitions_by_index.size(); n++ )
        {
            //            CellDefinition pCD = CellDefinitions_by_index[n]; 
            HypothesisRuleset pHRS = find_ruleset( pCD );

            String cell_type = pHRS.cell_type;
            try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
            {
                System.out.println( cell_type + " :: " );
                for( int k = 0; k < pHRS.rules.size(); k++ )
                {
                    String behavior = pHRS.rules.get( k ).behavior;
                    System.out.println( behavior + " : " );

                    double min_value = pHRS.rules.get( k ).min_value;
                    double max_value = pHRS.rules.get( k ).max_value;
                    double base_value = pHRS.rules.get( k ).base_value;
                    for( int i = 0; i < pHRS.rules.get( k ).signals.size(); i++ )
                    {
                        String signal = pHRS.rules.get( k ).signals.get( i );
                        System.out.println( signal + " " );
                        String response;
                        if( pHRS.rules.get( k ).responses.get( i ) == true )
                        {
                            response = "increases";
                        }
                        else
                        {
                            response = "decreases";
                        }
                        double half_max = pHRS.rules.get( k ).half_maxes.get( i );
                        double hill_power = pHRS.rules.get( k ).hill_powers.get( i );
                        boolean use_for_dead = false;

                        // output the rule 
                        bw.append( cell_type + "," + behavior + "," + min_value + "," + base_value + "," + max_value + "," + signal + ","
                                + response + "," + half_max + "," + hill_power + "," + use_for_dead + "\n" );
                    }
                    //                System.out.println( "\n"); 
                }
            }
        }

        /*
        0              1          2        3           4          5       6          7         8           9 
         Cell type, behavior, min value, base value, max value,   signal, direction, half-max, Hill power, dead
        */
        //        fs.close(); 

        System.out.println( "Done!" );
    }

    // v1 

    static void export_rules_csv_v1(String filename) throws Exception
    {
        File f = new File( filename );
        if( !f.exists() )
        {
            System.out.println( "Warning: Rules export file " + filename + " failed to open." );
            return;
        }
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            System.out.println( "Exporting rules to file " + filename + " (v1 format) ... " );

            for( CellDefinition pCD : CellDefinition.getCellDefinitions() )//int n=0; n < CellDefinitions_by_index.size(); n++ )
            {
                //            CellDefinition pCD = CellDefinitions_by_index[n]; 
                HypothesisRuleset pHRS = find_ruleset( pCD );

                String cell_type = pHRS.cell_type;

                for( int k = 0; k < pHRS.rules.size(); k++ )
                {
                    String behavior = pHRS.rules.get( k ).behavior;

                    double min_value = pHRS.rules.get( k ).min_value;
                    double max_value = pHRS.rules.get( k ).max_value;
                    double base_value = pHRS.rules.get( k ).base_value;
                    for( int i = 0; i < pHRS.rules.get( k ).signals.size(); i++ )
                    {
                        String signal = pHRS.rules.get( k ).signals.get( i );
                        String response;
                        double max_response = -9e99;
                        if( pHRS.rules.get( k ).responses.get( i ) == true )
                        {
                            response = "increases";
                            max_response = max_value;
                        }
                        else
                        {
                            response = "decreases";
                            max_response = min_value;
                        }
                        double half_max = pHRS.rules.get( k ).half_maxes.get( i );
                        double hill_power = pHRS.rules.get( k ).hill_powers.get( i );
                        boolean use_for_dead = false;

                        // output the rule 
                        bw.append( cell_type + "," + signal + "," + response + "," + behavior + "," // 0,1,2,3
                                + base_value + "," + max_response + "," + half_max + "," + hill_power + "," // 4,5, 6,7, 
                                + use_for_dead + "\n" );//std::endl; // 8 
                    }
                }
            }
        }

        /*
         Cell type, signal, direcxtion, behavior, base, max_response, half-max, hill , dead 
        */
        //        fs.close(); 

        System.out.println( "Done!" );
    }

    void export_rules_csv_v2(String filename) throws Exception
    {
        //        std::fstream fs( filename, std::ios::out );
        File f = new File( filename );
        if( !f.exists() )
        {
            System.out.println( "Warning: Rules export file " + filename + " failed to open." );
            return;
        }

        System.out.println( "Exporting rules to file " + filename + " (v2 format) ... " );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            for( CellDefinition pCD : CellDefinition.getCellDefinitions() )//int n=0; n < CellDefinitions_by_index.size(); n++ )
            {
                //            CellDefinition pCD = CellDefinitions_by_index[n]; 
                HypothesisRuleset pHRS = find_ruleset( pCD );

                String cell_type = pHRS.cell_type;

                for( int k = 0; k < pHRS.rules.size(); k++ )
                {
                    String behavior = pHRS.rules.get( k ).behavior;

                    double min_value = pHRS.rules.get( k ).min_value;
                    double max_value = pHRS.rules.get( k ).max_value;
                    double base_value = pHRS.rules.get( k ).base_value;
                    for( int i = 0; i < pHRS.rules.get( k ).signals.size(); i++ )
                    {
                        String signal = pHRS.rules.get( k ).signals.get( i );
                        String response;
                        double max_response = -9e99;
                        if( pHRS.rules.get( k ).responses.get( i ) == true )
                        {
                            response = "increases";
                            max_response = max_value;
                        }
                        else
                        {
                            response = "decreases";
                            max_response = min_value;
                        }
                        double half_max = pHRS.rules.get( k ).half_maxes.get( i );
                        double hill_power = pHRS.rules.get( k ).hill_powers.get( i );
                        boolean use_for_dead = false;

                        // output the rule 
                        bw.append( cell_type + "," + signal + "," + response + "," + behavior + "," // 0,1,2,3
                        // + base_value + "," 
                                + max_response + "," + half_max + "," + hill_power + "," // 4,5, 6,7, 
                                + use_for_dead + "\n" ); // 8 
                    }
                }
            }
        }

        /*
         Cell type, signal, direcxtion, behavior, base, max_response, half-max, hill , dead 
        */
        //        fs.close(); 
        System.out.println( "Done!" );
    }


    double[] UniformInUnitDisc()
    {
        double two_pi = 6.283185307179586;
        double theta = PhysiCellUtilities.UniformRandom(); // U(0,1)
        theta *= two_pi; // U(0,2*pi)
        double r = Math.sqrt( PhysiCellUtilities.UniformRandom() ); // sqrt( U(0,1) )
        return new double[] {r * Math.cos( theta ), r * Math.sin( theta ), 0.0};
    }

    double[] UniformInUnitSphere()
    {
        // reference: https://doi.org/10.1063/1.168311, adapting equation 13

        double two_pi = 6.283185307179586;

        double T = PhysiCellUtilities.UniformRandom();
        double sqrt_T = Math.sqrt( T );
        double sqrt_one_minus_T = 1.0;
        sqrt_one_minus_T -= T;
        sqrt_one_minus_T = Math.sqrt( sqrt_one_minus_T );

        double param1 = Math.pow( PhysiCellUtilities.UniformRandom(), 0.33333333333333333333333333333333333333 ); //  xi^(1/3), 
        double param2 = param1; // xi^(1/3)
        param2 *= 2.0; // 2 * xi^(1/3)
        param2 *= sqrt_T; // 2 * xi(1) * T^(1/2)
        param2 *= sqrt_one_minus_T; //  2 * xi(1) * T^(1/2) * (1-T)^(1/2)

        double theta = PhysiCellUtilities.UniformRandom(); // U(0,1)
        theta *= two_pi; // U(0,2*pi)

        return new double[] {param2 * Math.sin( theta ), param2 * Math.cos( theta ), param1 * ( 1 - 2 * T )};
    }

    double[] UniformInAnnulus(double r1, double r2)
    {
        double two_pi = 6.283185307179586;

        double theta = PhysiCellUtilities.UniformRandom();
        theta *= two_pi;
        double r1_2 = r1 * r1;
        double r2_2 = r2 * r2;

        double r = Math.sqrt( r1_2 + ( r2_2 - r1_2 ) * PhysiCellUtilities.UniformRandom() );
        double x = r * Math.cos( theta );
        double y = r * Math.sin( theta );
        return new double[] {x, y, 0.0};
    }

    double[] UniformInShell(double r1, double r2)
    {
        double two_pi = 6.283185307179586;

        double T = PhysiCellUtilities.UniformRandom();
        double sqrt_T = Math.sqrt( T );
        double sqrt_one_minus_T = 1.0;
        sqrt_one_minus_T -= T;
        sqrt_one_minus_T = Math.sqrt( sqrt_one_minus_T );

        double param1 = Math.pow( PhysiCellUtilities.UniformRandom(), 0.33333333333333333333333333333333333333 ); //  xi^(1/3), 
        // param1 *= (r2-r1); 
        // param1 += r1; 
        double param2 = param1; // xi^(1/3)
        param2 *= 2.0; // 2 * xi^(1/3)
        param2 *= sqrt_T; // 2 * xi(1) * T^(1/2)
        param2 *= sqrt_one_minus_T; //  2 * xi(1) * T^(1/2) * (1-T)^(1/2)

        double theta = PhysiCellUtilities.UniformRandom(); // U(0,1)
        theta *= two_pi; // U(0,2*pi)

        return new double[] {param2 * Math.sin( theta ), param2 * Math.cos( theta ), param1 * ( 1 - 2 * T )};
    }

    public static void setup_cell_rules(Model model )
    {
        // setup 
        intialize_hypothesis_rulesets(); 

        // load rules 
        //        parse_rules_from_pugixml(); 

        // display rules to screen
        display_hypothesis_rulesets( );

        // save annotations 
//        save_annotated_detailed_English_rules(); 
//        save_annotated_detailed_English_rules_HTML(); 
//        save_annotated_English_rules(); 
//        save_annotated_English_rules_HTML(); 

        // save dictionaries 
        //        String dictionary_file = "./" + PhysiCellSettings.folder + "/dictionaries.txt";
        //        std::ofstream dict_of( dictionary_file , std::ios::out ); 

        //        display_signal_dictionary( dict_of ); // done 
        //        display_behavior_dictionary( dict_of ); // done 
        //        dict_of.close(); 

        // save rules (v1)
        //        String rules_file = PhysiCellSettings.folder + "/cell_rules.csv"; 
        //        export_rules_csv_v1( rusles_file ); 
    }

}
