package ru.biosoft.physicell.core;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
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
    public Map<CellDefinition, HypothesisRuleset> hypothesisRulesets = new HashMap<>();

    void add_hypothesis_ruleset(Model model, CellDefinition cd)
    {
        HypothesisRuleset ruleset = new HypothesisRuleset();
        ruleset.sync( model, cd );
        hypothesisRulesets.put( cd, ruleset );
    }

    public void initRulesets(Model model)
    {
        hypothesisRulesets.clear();

        for( CellDefinition cd : model.getCellDefinitions() )
        {
            add_hypothesis_ruleset( model, cd );
        }
    }

    HypothesisRuleset access_ruleset(CellDefinition cd)
    {
        return hypothesisRulesets.get( cd );
    }

    public HypothesisRuleset findRuleset(CellDefinition cd)
    {
        return ( hypothesisRulesets.get( cd ) );
    }

    public String display_hypothesis_rulesets(Model model)
    {
        StringBuilder sb = new StringBuilder();
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            sb.append( hypothesisRulesets.get( cd ).display() );
        }
        return sb.toString();
    }

    String detailed_display_hypothesis_rulesets(Model model)
    {
        StringBuilder sb = new StringBuilder();
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            sb.append( hypothesisRulesets.get( cd ).detailedDisplay() );
        }
        return sb.toString();
    }

    public static void addRule(Model model, String type, String signal, String behavior, String response) throws Exception
    {
        CellDefinition cd = model.getCellDefinition( type );
        if( cd == null )
            System.out.println( "Warning: Attempted to add rule for " + type + ", but no cell definition found for this type." );

        HypothesisRuleset pHRS = model.rules.findRuleset( cd );
        if( pHRS == null )
            System.out.println( "Warning: Attempted to add rule for " + type + ", but no hypothesis ruleset found for this type." );

        if( pHRS.findBehavior( behavior ) != null )
        {
            if( !pHRS.get( behavior ).behavior.equals( behavior ) )
            {
                pHRS.get( behavior ).behavior = behavior;
                System.out.println( "wha?" );
            }
        }
        pHRS.addBehavior( model, behavior );
        pHRS.get( behavior ).addSignal( signal, response );
    }

    public static void addRule(Model model, String type, String signal, String behavior, String response, boolean useForDead)
            throws Exception
    {
        CellDefinition cd = model.getCellDefinition( type );
        if( cd == null )
            throw new Exception( "Warning: Attempted to add rule for " + type + ", but no cell definition found for this type." );

        HypothesisRuleset pHRS = model.rules.findRuleset( cd );
        if( pHRS == null )
            throw new Exception( "Warning: Attempted to add rule for " + type + ", but no hypothesis ruleset found for this type." );

        if( pHRS.findBehavior( behavior ) != null )
        {
            if( !pHRS.get( behavior ).behavior.equals( behavior ) )
            {
                pHRS.get( behavior ).behavior = behavior;
                System.out.println( "wha?" );
            }
        }
        pHRS.addBehavior( model, behavior );
        pHRS.get( behavior ).addSignal( signal, response );
        // set dead flag
        int n = pHRS.get( behavior ).findSignal( signal );
        pHRS.get( behavior ).appliesToDead.set( n, useForDead );
    }

    void set_hypothesis_parameters(Model model, String type, String signal, String behavior, double max, double hillPower)
            throws Exception
    {
        CellDefinition cd = model.getCellDefinition( type );
        if( cd == null )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + type
                    + ", but no cell definition found for this type." );

        if( findRuleset( cd ) == null )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + type
                    + ", but no behavior ruleset for this cell type." );

        if( hypothesisRulesets.get( cd ).findBehavior( behavior ) == null )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + type
                    + ", but the cell type has no rules for this behavior." );
        // System.out.println( "Error. No " + behavior + " rules for " + cell_type + ". Ignoring."+ std::endl; return;

        if( hypothesisRulesets.get( cd ).get( behavior ).findSignal( signal ) < 0 )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " modulated by " + signal + " in " + type
                    + ", but the cell type's behavior does not vary with this signal." );

        hypothesisRulesets.get( cd ).get( behavior ).setHalfMax( signal, max );
        hypothesisRulesets.get( cd ).get( behavior ).setHillPower( signal, hillPower );
    }

    public static void setBehaviorParameters(Model model, String type, String behavior, double min, double max) throws Exception
    {
        CellDefinition cd = model.getCellDefinition( type );
        if( cd == null )
            throw new Exception(
                    "Warning: Attempted to set parameters for " + behavior + " in " + type + ", but the cell definition is not found." );

        if( model.rules.findRuleset( cd ) == null )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + type
                    + ", but there is no hypothesis ruleset for this cell type." );

        if( model.rules.hypothesisRulesets.get( cd ).findBehavior( behavior ) == null )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + type
                    + ", but there is no rules for this behavior for this cell type." );

        model.rules.hypothesisRulesets.get( cd ).get( behavior ).minValue = min;
        model.rules.hypothesisRulesets.get( cd ).get( behavior ).maxValue = max;
    }

    public static void setBehaviorParameters(Model model, String type, String behavior, double min, double base, double max)
            throws Exception
    {
        CellDefinition pCD = model.getCellDefinition( type );
        if( pCD == null )
            throw new Exception(
                    "Warning: Attempted to set parameters for " + behavior + " in " + type + ", but the cell definition is not found." );

        if( model.rules.findRuleset( pCD ) == null )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + type
                    + ", but there is no hypothesis ruleset for this cell type." );

        if( model.rules.hypothesisRulesets.get( pCD ).findBehavior( behavior ) == null )
            throw new Exception( "Warning: Attempted to set parameters for " + behavior + " in " + type
                    + ", but there is no rules for this behavior for this cell type." );

        model.rules.hypothesisRulesets.get( pCD ).get( behavior ).minValue = min;
        model.rules.hypothesisRulesets.get( pCD ).get( behavior ).maxValue = max;
        model.rules.hypothesisRulesets.get( pCD ).get( behavior ).baseValue = base;
    }


    public static void setBehaviorBaseValue(Model model, String type, String behavior, double base) throws Exception
    {
        CellDefinition cd = model.getCellDefinition( type );
        if( cd == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but the cell definition is not found." );

        if( model.rules.findRuleset( cd ) == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but no hypothesis ruleset found for this cell type." );

        if( model.rules.hypothesisRulesets.get( cd ).findBehavior( behavior ) == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but no rule for this behavior found for this cell type." );

        model.rules.hypothesisRulesets.get( cd ).get( behavior ).baseValue = base;
    }

    public static void setBehaviorMinValue(Model model, String type, String behavior, double value) throws Exception
    {
        CellDefinition cd = model.getCellDefinition( type );
        if( cd == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but the cell definition is not found." );

        if( model.rules.findRuleset( cd ) == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but no hypothesis ruleset found for this cell type." );

        if( model.rules.hypothesisRulesets.get( cd ).findBehavior( behavior ) == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but no rule for this behavior found for this cell type." );

        model.rules.hypothesisRulesets.get( cd ).get( behavior ).minValue = value;
    }

    public static void set_behavior_max_value(Model model, String type, String behavior, double max) throws Exception
    {
        CellDefinition pCD = model.getCellDefinition( type );
        if( pCD == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but the cell definition is not found." );

        if( model.rules.findRuleset( pCD ) == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but no hypothesis ruleset found for this cell type." );

        if( model.rules.hypothesisRulesets.get( pCD ).findBehavior( behavior ) == null )
            throw new Exception( "Warning: Attempted to set base parameter for " + behavior + " in " + type
                    + ", but no rule for this behavior found for this cell type." );

        model.rules.hypothesisRulesets.get( pCD ).get( behavior ).maxValue = max;
    }

    void applyRuleset(Cell cell) throws Exception
    {
        CellDefinition cd = cell.getModel().getCellDefinition( cell.typeName );
        HypothesisRuleset set = hypothesisRulesets.get( cd );
        if( set != null )
            set.apply( cell );
    }

    void rule_phenotype_function(Cell cell, Phenotype phenotype, double dt) throws Exception
    {
        applyRuleset( cell );

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

    static String csv_strings_to_English_v2(String[] strings, boolean include_cell_header)
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

    void parseCSVRule0(Model model, String[] input) throws Exception
    {
        boolean skip = input.length != 10;

        for( int n = 0; n < input.length; n++ )
        {
            if( input[n].isEmpty() )
                skip = true;
        }
        if( skip )
        {
            System.out.println( "Warning: Misformed rule (likely from an empty rules file). Skipping." );
            return;
        }

        String temp = csv_strings_to_English( input, false );

        // string portions of the rule
        String type = input[0];
        String behavior = input[1];
        String signal = input[5];
        String response = input[6];

        // numeric portions of the rule 
        double min = Double.parseDouble( input[2] );
        double base = Double.parseDouble( input[3] );
        double max = Double.parseDouble( input[4] );

        double halfMax = Double.parseDouble( input[7] );
        double hillPower = Double.parseDouble( input[8] );
        boolean useForDead = parseBoolean( input[9] );

        //        System.out.println( "Adding rule for " + type + " cells:\n\t" );
        //        System.out.println( temp );

        // add_rule(cell_type,signal,behavior,response);  
        addRule( model, type, signal, behavior, response, useForDead );

        set_hypothesis_parameters( model, type, signal, behavior, halfMax, hillPower );

        // compare to base behavior value in cell def for discrepancies 
        CellDefinition cd = model.getCellDefinition( type );
        double refBaseValue = model.getSignals().getSingleBaseBehavior( model, cd, behavior );
        if( Math.abs( refBaseValue - base ) > 1e-15 )
        {
            throw new Exception( "Error: Base value for " + behavior + " in cell type " + type + "\n" + "       has base value " + base
                    + " in the rule, " + "\n" + "but base value " + refBaseValue + " in the cell definition." + "\n"
                    + "       Fix this discrepancy to continue the model." );
        }

        // get_single_base_behavior

        setBehaviorParameters( model, type, behavior, min, base, max );

        /*
        // set "applies to dead"
        
        int n = access_ruleset(pCD)[behavior].find_signal(signal); 
        access_ruleset(pCD)[behavior].applies_to_dead_cells[n] = use_for_dead; 
        */
    }

    void parseCSVRule0(Model model, String input) throws Exception
    {
        String[] tokenized = input.split( "," );
        if( tokenized.length != 10 )
            tokenized = input.split( "\t" );
        parseCSVRule0( model, tokenized );
    }

    void parseCSVRules0(Model model, InputStream stream) throws Exception
    {
        try (BufferedReader br = new BufferedReader( new InputStreamReader( stream ) ))
        {
            String line = br.readLine();
            while( line != null )
            {
                if( line.length() > 0 )
                    parseCSVRule0( model, line );
            }
            line = br.readLine();
        }
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

    void parseCSVRule1(Model model, String[] input) throws Exception
    {
        boolean skip = false;
        if( input.length != 9 )
            skip = true;

        for( int n = 0; n < input.length; n++ )
        {
            if( input[n].isEmpty() )
                skip = true;
        }
        if( skip )
        {
            System.out.println( "Warning: Misformed rule (likely from an empty rules file). Skipping." );
            for( int n = 0; n < input.length; n++ )
                System.out.println( n + " : " + input[n] );
            return;
        }

        String temp = csv_strings_to_English_v1( input, false ); // need a v1 version of this

        // string portions of the rule
        String type = input[0];
        String signal = input[1];
        String response = input[2];
        String behavior = input[3];

        // numeric portions of the rule 
        // double min_value  = std::atof( input[2].c_str() );

        double base = Double.parseDouble( input[4] );
        double maxResponse = Double.parseDouble( input[5] );

        // hmm from here 
        // double max_value  = std::atof( input[4].c_str() ); 

        double halfMax = Double.parseDouble( input[6] );
        double hillPower = Double.parseDouble( input[7] );
        boolean useForDead = parseBoolean( input[8] );

        //        System.out.println( "Adding rule for " + type + " cells:\n\t" );
        //        System.out.println( temp );

        // add_rule(cell_type,signal,behavior,response);  
        addRule( model, type, signal, behavior, response, useForDead );
        set_hypothesis_parameters( model, type, signal, behavior, halfMax, hillPower );

        // compare to base behavior value in cell def for discrepancies 
        CellDefinition cd = model.getCellDefinition( type );
        double refBaseValue = model.getSignals().getSingleBaseBehavior( model, cd, behavior );
        if( Math.abs( refBaseValue - base ) > 1e-15 )
        {
            throw new Exception( "Error: Base value for " + behavior + " in cell type " + type + "\n" + "       has base value " + base
                    + " in the rule, " + "\n"// + std::endl 
                    + "but base value " + refBaseValue + " in the cell definition." + "\n"
                    + "       Fix this discrepancy to continue the model." );
        }

        setBehaviorBaseValue( model, type, behavior, base );

        if( response.equals( "increases" ) )
        {
            set_behavior_max_value( model, type, behavior, maxResponse );
        }
        else
        {
            setBehaviorMinValue( model, type, behavior, maxResponse );
        }
    }

    void parseCSVRule1(Model model, String input) throws Exception
    {
        String[] tokenized_string = input.split( "," );
        if( tokenized_string.length != 9 )
            tokenized_string = input.split( ",\t" );

        // check for comment 
        if( tokenized_string[0].charAt( 0 ) == '/' && tokenized_string[0].charAt( 1 ) == '/' )
        {
            System.out.println( "Skipping commented rule (" + input + ")" );
            return;
        }
        parseCSVRule1( model, tokenized_string );
    }

    void parseCSVRules1(Model model, String filename) throws Exception
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
                    parseCSVRule1( model, line );
                line = br.readLine();
            }
        }
    }


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

    static void parseCSVRule2(Model model, String[] input) throws Exception
    {
        Rules rules = model.rules;
        boolean skip = false;
        if( input.length != 8 )
            skip = true;

        for( int n = 0; n < input.length; n++ )
        {
            if( input[n].isEmpty() )
                skip = true;
        }
        if( skip )
        {
            System.out.println( "Warning: Misformed rule (likely from an empty rules file). Skipping." );
            for( int n = 0; n < input.length; n++ )
                System.out.println( n + " : " + input[n] );
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
        double max_response = Double.parseDouble( input[4] );

        // hmm from here 
        // double max_value  = std::atof( input[4].c_str() ); 

        double half_max = Double.parseDouble( input[5] );
        double hill_power = Double.parseDouble( input[6] );
        boolean use_for_dead = parseBoolean( input[7] );

        //        System.out.println( "Adding rule for " + cell_type + " cells:\n\t" );
        //        System.out.println( temp );

        // add_rule(cell_type,signal,behavior,response);  
        addRule( model, cell_type, signal, behavior, response, use_for_dead );
        rules.set_hypothesis_parameters( model, cell_type, signal, behavior, half_max, hill_power );

        // compare to base behavior value in cell def for discrepancies 
        CellDefinition pCD = model.getCellDefinition( cell_type );
        double ref_base_value = model.getSignals().getSingleBaseBehavior( model, pCD, behavior );
        rules.setBehaviorBaseValue( model, cell_type, behavior, ref_base_value );
        if( response.equals( "increases" ) )
            rules.set_behavior_max_value( model, cell_type, behavior, max_response );
        else
            rules.setBehaviorMinValue( model, cell_type, behavior, max_response );
    }

    public static boolean parseBoolean(String val)
    {
        return Double.parseDouble( val ) > 0;
    }

    static void parseCSVRule2(Model model, String input) throws Exception
    {
        if( input.startsWith( "//" ) )
            return;
        String[] tokenized_string = input.split( "," );
        if( tokenized_string.length != 8 )
            tokenized_string = input.split( "\t" );
        //        if( tokenized_string[0].charAt( 0 ) == '/' && tokenized_string[0].charAt( 1 ) == '/' )
        //            throw new Exception( "Skipping commented rule (" + input + ")" );
        parseCSVRule2( model, tokenized_string );
    }

    public static void parseCSVRules2(Model model, String filePath) throws Exception
    {
        parseCSVRules2( model, new FileInputStream( new File( filePath ) ) );
    }

    public static void parseCSVRules2(Model model, InputStream stream) throws Exception
    {
        try (BufferedReader br = new BufferedReader( new InputStreamReader( stream ) ))
        {
            String line = br.readLine();
            while( line != null )
            {
                if( line.length() > 0 )
                    parseCSVRule2( model, line );
                line = br.readLine();
            }
        }
    }

    public static String stream_annotated_English_rules(Model model)
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Cell Hypothesis Rules\n\n" );
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            HypothesisRuleset pHRS = model.rules.findRuleset( cd );
            sb.append( "In " + pHRS.type + " cells:\n" );
            for( int k = 0; k < pHRS.rules.size(); k++ )
                sb.append( pHRS.rules.get( k ).English_display() );
            sb.append( "\n" );
        }
        return sb.toString();
    }

    public static String stream_annotated_English_rules_HTML(Model model)
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "<html>\n<body><h1>Cell Hypothesis Rules</h1>\n" );
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            sb.append( "<p>" );
            HypothesisRuleset pHRS = model.rules.findRuleset( cd );
            sb.append( "In " + pHRS.type + " cells:\n" );
            sb.append( "<ul>\n" );
            for( int k = 0; k < pHRS.rules.size(); k++ )
                sb.append( pHRS.rules.get( k ).English_display_HTML() );
            sb.append( "</ul>\n</p>\n" );
        }
        sb.append( "</body>\n</html>\n" );
        return sb.toString();
    }

    static void save_annotated_English_rules(Model model, String filename)
    {
        //        String filename = PhysiCellSettings.folder + "/rules.txt";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_English_rules( model );
        //        of.close(); 
    }

    public static String stream_annotated_detailed_English_rules(Model model)
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Cell Hypothesis Rules (detailed)\n\n" );
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            HypothesisRuleset pHRS = model.rules.findRuleset( cd );
            sb.append( "In " + pHRS.type + " cells:\n" );
            for( int k = 0; k < pHRS.rules.size(); k++ )
            {
                pHRS.rules.get( k ).English_detailed_display();
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }

    static String stream_annotated_detailed_English_rules_HTML(Model model)
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "<html>\n<body><h1>Cell Hypothesis Rules (detailed)</h1>\n" );
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            sb.append( "<p>" );
            HypothesisRuleset pHRS = model.rules.findRuleset( cd );
            sb.append( "In " + pHRS.type + " cells:\n" );
            sb.append( "<ul>" );// + std::endl; 
            for( int k = 0; k < pHRS.rules.size(); k++ )
            {
                sb.append( pHRS.rules.get( k ).English_detailed_display_HTML() );
            }
            sb.append( "</ul>\n</p>\n" );
        }
        sb.append( "</body>\n</html>\n" );
        return sb.toString();
    }

    static void save_annotated_detailed_English_rules(Model model, String filename) throws Exception
    {
        //        String filename = PhysiCell_settings.folder + "/detailed_rules.txt";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_detailed_English_rules( model );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            bw.write( text );
        }
        //        of.close(); 
    }

    static void save_annotated_detailed_English_rules_HTML(Model model, String filename) throws Exception
    {
        //        String filename = PhysiCell_settings.folder + "/detailed_rules.html";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_detailed_English_rules_HTML( model );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            bw.write( text );
        }
        //        of.close(); 
    }

    static void save_annotated_English_rules_HTML(Model model, String filename) throws Exception
    {
        //        String filename = PhysiCell_settings.folder + "/rules.html";
        //        std::ofstream of( filename , std::ios::out );
        String text = stream_annotated_English_rules_HTML( model );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
        {
            bw.write( text );
        }
        //        of.close(); 
    }

    // v0 version 

    static void export_rules_csv_v0(Model model, String filename) throws Exception
    {
        File f = new File( filename );
        if( !f.exists() )
        {
            System.out.println( "Warning: Rules export file " + filename + " failed to open." );
            return;
        }

        System.out.println( "Exporting rules to file " + filename + " (v0 format) ... " );

        for( CellDefinition pCD : model.getCellDefinitions() )// n=0; n < CellDefinitions_by_index.size(); n++ )
        {
            //            CellDefinition pCD = CellDefinitions_by_index[n]; 
            HypothesisRuleset pHRS = model.rules.findRuleset( pCD );

            String cell_type = pHRS.type;
            try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( filename ) ) ))
            {
                System.out.println( cell_type + " :: " );
                for( int k = 0; k < pHRS.rules.size(); k++ )
                {
                    String behavior = pHRS.rules.get( k ).behavior;
                    System.out.println( behavior + " : " );

                    double min_value = pHRS.rules.get( k ).minValue;
                    double max_value = pHRS.rules.get( k ).maxValue;
                    double base_value = pHRS.rules.get( k ).baseValue;
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
                        double half_max = pHRS.rules.get( k ).halfMaxes.get( i );
                        double hill_power = pHRS.rules.get( k ).hillPowers.get( i );
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

    static void export_rules_csv_v1(Model model, String filename) throws Exception
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

            for( CellDefinition pCD : model.getCellDefinitions() )//int n=0; n < CellDefinitions_by_index.size(); n++ )
            {
                //            CellDefinition pCD = CellDefinitions_by_index[n]; 
                HypothesisRuleset pHRS = model.rules.findRuleset( pCD );

                String cell_type = pHRS.type;

                for( int k = 0; k < pHRS.rules.size(); k++ )
                {
                    String behavior = pHRS.rules.get( k ).behavior;

                    double min_value = pHRS.rules.get( k ).minValue;
                    double max_value = pHRS.rules.get( k ).maxValue;
                    double base_value = pHRS.rules.get( k ).baseValue;
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
                        double half_max = pHRS.rules.get( k ).halfMaxes.get( i );
                        double hill_power = pHRS.rules.get( k ).hillPowers.get( i );
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

    void export_rules_csv_v2(String filename, Model model) throws Exception
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
            for( CellDefinition pCD : model.getCellDefinitions() )//int n=0; n < CellDefinitions_by_index.size(); n++ )
            {
                //            CellDefinition pCD = CellDefinitions_by_index[n]; 
                HypothesisRuleset pHRS = findRuleset( pCD );

                String cell_type = pHRS.type;

                for( int k = 0; k < pHRS.rules.size(); k++ )
                {
                    String behavior = pHRS.rules.get( k ).behavior;

                    double min_value = pHRS.rules.get( k ).minValue;
                    double max_value = pHRS.rules.get( k ).maxValue;
                    double base_value = pHRS.rules.get( k ).baseValue;
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
                        double half_max = pHRS.rules.get( k ).halfMaxes.get( i );
                        double hill_power = pHRS.rules.get( k ).hillPowers.get( i );
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


    double[] UniformInUnitDisc(Model model)
    {
        double two_pi = 6.283185307179586;
        double theta = model.rng.UniformRandom(); // U(0,1)
        theta *= two_pi; // U(0,2*pi)
        double r = Math.sqrt( model.rng.UniformRandom() ); // sqrt( U(0,1) )
        return new double[] {r * Math.cos( theta ), r * Math.sin( theta ), 0.0};
    }

    double[] UniformInUnitSphere(Model model)
    {
        // reference: https://doi.org/10.1063/1.168311, adapting equation 13

        double two_pi = 6.283185307179586;

        double T = model.rng.UniformRandom();
        double sqrt_T = Math.sqrt( T );
        double sqrt_one_minus_T = 1.0;
        sqrt_one_minus_T -= T;
        sqrt_one_minus_T = Math.sqrt( sqrt_one_minus_T );

        double param1 = Math.pow( model.rng.UniformRandom(), 0.33333333333333333333333333333333333333 ); //  xi^(1/3), 
        double param2 = param1; // xi^(1/3)
        param2 *= 2.0; // 2 * xi^(1/3)
        param2 *= sqrt_T; // 2 * xi(1) * T^(1/2)
        param2 *= sqrt_one_minus_T; //  2 * xi(1) * T^(1/2) * (1-T)^(1/2)

        double theta = model.rng.UniformRandom(); // U(0,1)
        theta *= two_pi; // U(0,2*pi)

        return new double[] {param2 * Math.sin( theta ), param2 * Math.cos( theta ), param1 * ( 1 - 2 * T )};
    }

    double[] UniformInAnnulus(Model model, double r1, double r2)
    {
        double two_pi = 6.283185307179586;

        double theta = model.rng.UniformRandom();
        theta *= two_pi;
        double r1_2 = r1 * r1;
        double r2_2 = r2 * r2;

        double r = Math.sqrt( r1_2 + ( r2_2 - r1_2 ) * model.rng.UniformRandom() );
        double x = r * Math.cos( theta );
        double y = r * Math.sin( theta );
        return new double[] {x, y, 0.0};
    }

    double[] UniformInShell(Model model, double r1, double r2)
    {
        double two_pi = 6.283185307179586;

        double T = model.rng.UniformRandom();
        double sqrt_T = Math.sqrt( T );
        double sqrt_one_minus_T = 1.0;
        sqrt_one_minus_T -= T;
        sqrt_one_minus_T = Math.sqrt( sqrt_one_minus_T );

        double param1 = Math.pow( model.rng.UniformRandom(), 0.33333333333333333333333333333333333333 ); //  xi^(1/3), 
        // param1 *= (r2-r1); 
        // param1 += r1; 
        double param2 = param1; // xi^(1/3)
        param2 *= 2.0; // 2 * xi^(1/3)
        param2 *= sqrt_T; // 2 * xi(1) * T^(1/2)
        param2 *= sqrt_one_minus_T; //  2 * xi(1) * T^(1/2) * (1-T)^(1/2)

        double theta = model.rng.UniformRandom(); // U(0,1)
        theta *= two_pi; // U(0,2*pi)

        return new double[] {param2 * Math.sin( theta ), param2 * Math.cos( theta ), param1 * ( 1 - 2 * T )};
    }

    public static void setupRules(Model model)
    {
        model.rules.initRulesets( model );
    }
}
