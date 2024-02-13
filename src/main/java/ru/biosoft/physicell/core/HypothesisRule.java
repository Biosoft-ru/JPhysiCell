package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.List;
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
public class HypothesisRule
{
    Map<String, Integer> signals_map;

    String cell_type;
    CellDefinition pCell_Definition;

    String behavior;
    double base_value;
    double max_value;
    double min_value;

    List<String> signals;
    List<Boolean> responses;
    List<Double> half_maxes;
    List<Double> hill_powers;
    List<Boolean> applies_to_dead_cells;

    List<String> up_signals;
    List<Double> up_half_maxes;
    List<Double> up_hill_powers;
    List<Boolean> up_applies_to_dead_cells;

    List<String> down_signals;
    List<Double> down_half_maxes;
    List<Double> down_hill_powers;
    List<Boolean> down_applies_to_dead_cells;

    HypothesisRule()
    {
        signals_map.clear();

        behavior = "none";
        base_value = 1.0;
        max_value = 10.0;
        min_value = 0.1;

        signals = new ArrayList<>();//String[0];//.resize( 0 );
        responses = new ArrayList<>();//boolean[0];
        half_maxes = new ArrayList<>();//new double[0];
        hill_powers = new ArrayList<>();//new double[0];
        applies_to_dead_cells = new ArrayList<>();//new boolean[0];

        up_signals = new ArrayList<>();//= new String[0];
        up_half_maxes = new ArrayList<>();//= new double[0];
        up_hill_powers = new ArrayList<>();//= new double[0];
        up_applies_to_dead_cells = new ArrayList<>();//= new boolean[0];

        down_signals = new ArrayList<>();//= new String[0];
        down_half_maxes = new ArrayList<>();//= new double[0];
        down_hill_powers = new ArrayList<>();//= new double[0];
        down_applies_to_dead_cells = new ArrayList<>();//= new boolean[0];

        cell_type = "none";
        pCell_Definition = null;
    }

    String convert_bool_to_response(boolean input)
    {
        if( input )
        {
            return "increases";
        }
        return "decreases";
    }

    /*
    double multivariate_Hill_response_function( double[] signals, double[] half_maxes , double[] hill_powers )
    {
        double temp1 = 0.0; 
        double temp2 = 0.0; 
        double temp3 = 0.0; 
        // create the generalized (s^h), stored in temp1; 
        for( int j=0 ; j < signals.size(); j++ )
        {
            temp2 = signals[j];     // s
            temp2 /= half_maxes[j]; // s/s_half 
            temp3 = pow( temp2 , hill_powers[j] ); // (s/s_half)^h 
            temp1 += temp3; 
        }
        temp2 = temp1;   // numerator (S^h)
        temp1 += 1.0;    // denominator (1+S^h)
        temp2 /= temp1;  // numerator/denominator = S^h / (1+S^h)
        return temp2; 
    }
    
    double multivariate_linear_response_function( double[] signals, double[] min_thresholds , double[] max_thresholds ) 
    {
        double output = 0.0; 
    
        for( int j=0 ; j < signals.length; j++ )
        { output += linear_response_function( signals[j] , min_thresholds[j], max_thresholds[j] ); }
    
        if( output > 1.0 )
        { return 1.0; }
    
        return output; 
    }
    */

    String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "For cell type " + cell_type + ": \n" );
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( signals.get( j ) + " " + convert_bool_to_response( responses.get( j ) ) + " " + behavior + "\n" );
        }
        return sb.toString();
    }

    String reduced_display()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( signals.get( j ) + " " + convert_bool_to_response( responses.get( j ) ) + " " + behavior + "\n" );
        }

        return sb.toString();
    }

    String detailed_display()
    {
        StringBuilder sb = new StringBuilder();
        // sb.append( "For cell type " + cell_type + ": " + std::endl; 
        sb.append( behavior + " is modulated from " + min_value + " to " + max_value + " with a base value of " + base_value + "\n" );
        sb.append( "--------------------------------------------------------\n" );
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( "\t" + signals.get( j ) + " " + convert_bool_to_response( responses.get( j ) ) + " " + behavior + " with half-max "
                    + half_maxes.get( j ) + " and Hill power " + hill_powers.get( j ) + "." );
            if( applies_to_dead_cells.get( j ) == true )
            {
                sb.append( " Rule applies to dead cells." );
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }

    String English_detailed_display()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( signals.get( j ) + " " );
            if( responses.get( j ) == true )
            {
                sb.append( "increases " );
            }
            else
            {
                sb.append( "decreases " );
            }
            sb.append( behavior + " from " + base_value + " towards " );
            if( responses.get( j ) == true )
            {
                sb.append( max_value );
            }
            else
            {
                sb.append( min_value );
            }
            sb.append( " with a Hill response, with half-max " + half_maxes.get( j ) );
            sb.append( " and Hill power " + hill_powers.get( j ) + "." );
            if( applies_to_dead_cells.get( j ) == true )
            {
                sb.append( " Rule applies to dead cells." );
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }

    String English_detailed_display_HTML()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( "<li>" + signals.get( j ) + " " );
            if( responses.get( j ) == true )
            {
                sb.append( "increases " );
            }
            else
            {
                sb.append( "decreases " );
            }
            sb.append( behavior + " from " + base_value + " towards " );
            if( responses.get( j ) == true )
            {
                sb.append( max_value );
            }
            else
            {
                sb.append( min_value );
            }
            sb.append( " with a Hill response, with half-max " + half_maxes.get( j ) );
            sb.append( " and Hill power " + hill_powers.get( j ) + "." );
            if( applies_to_dead_cells.get( j ) == true )
            {
                sb.append( " Rule applies to dead cells." );
            }
            sb.append( "</li>\n" );
        }
        return sb.toString();
    }

    String English_display()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( signals.get( j ) + " " );
            if( responses.get( j ) == true )
            {
                sb.append( "increases " );
            }
            else
            {
                sb.append( "decreases " );
            }
            sb.append( behavior + "\n" );
        }
        return sb.toString();
    }

    String English_display_HTML()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( "<li>" + signals.get( j ) + " " );
            if( responses.get( j ) == true )
            {
                sb.append( "increases " );
            }
            else
            {
                sb.append( "decreases " );
            }
            sb.append( behavior + "</li>\n" );
        }
        return sb.toString();
    }

    void add_signal(String signal, double half_max, double hill_power, String response) throws Exception
    {
        // check: is this a valid signal? (is it in the dictionary?)
        if( SignalBehavior.find_signal_index( signal ) < 0 )
        {
            throw new Exception( "Warning! Attempted to add signal " + signal + " which is not in the dictionary."
                    + "Either fix your model or add the missing signal to the simulation." );
        }

        // check to see if it's already there 
        int n = find_signal( signal );

        // if so, then just warn and exit.  
        if( n > -1 )
        {
            System.out.println( "Warning! Signal " + signal + " was already part of the rule. Ignoring input." );
            return;
        }

        // add the signal; 
        signals_map.put( signal, signals_map.size() );

        boolean bResponse = false; // true if up-regulate, false if down
        if( response == "increase" || response == "increases" || response == "promotes" )
        {
            bResponse = true;
        }

        signals.add( signal );
        half_maxes.add( half_max );
        hill_powers.add( hill_power );
        responses.add( bResponse );
        applies_to_dead_cells.add( false );

        // separate into up and down for our convenience 
        if( bResponse == true )
        {
            up_signals.add( signal );
            up_half_maxes.add( half_max );
            up_hill_powers.add( hill_power );
            up_applies_to_dead_cells.add( false );
        }
        else
        {
            down_signals.add( signal );
            down_half_maxes.add( half_max );
            down_hill_powers.add( hill_power );
            down_applies_to_dead_cells.add( false );
        }
    }

    void add_signal(String signal, String response) throws Exception
    {
        add_signal( signal, 0.5, 3.0, response );
    }

    double evaluate( double[] signal_values , boolean dead )
    {
    // create signals 
    List<Double> up_signal = new ArrayList<>();//(0,0); //  up_signals.size() , 0.0 ); 
    List<Double> down_signal = new ArrayList<>();//(0,0); //  down_signals.size() , 0.0 ); 

        boolean apply_rule = false; 

        // need to modify to evaluate if cell is live or if the rule is allowed for dead cells 

        for( int j = 0; j < signal_values.length; j++ )
        { 
            if( applies_to_dead_cells.get(j) == true || dead == false )
            {
                if( responses.get(j) )
                {
                    up_signal.add( signal_values[j] );
                }
                else
                {
                    down_signal.add( signal_values[j] );
                }

                apply_rule = true; 
            }
            else
            {
                // new oin sep 7 , 2022
                if( responses.get(j) )
                {
                    up_signal.add( 0.0 );
                }
                else
                {
                    down_signal.add( 0.0 );
                }
            }
        }

        // March 27, 2023
        // if none of the rules apply, the return an absurdly low value 
        // to signal that the parameter value shoudl not be written 
        if( apply_rule == false )
        { return -9e99; }

        // up-regulation part 
        double HU = BasicSignaling.multivariate_Hill_response_function( up_signal, up_half_maxes, up_hill_powers );
        double U = base_value + (max_value-base_value)*HU; 

        // then the down-regulation part 
        double DU = BasicSignaling.multivariate_Hill_response_function( down_signal, down_half_maxes, down_hill_powers );
        double output = U + (min_value-U)*DU; 

        return output; 
    }

    double evaluate(double[] signal_values)
    {
        return evaluate( signal_values, true );
    }

    double evaluate(Cell pCell) throws Exception
    {
    // construct signal vector 
    double[] signal_values = new double[signals.size()];//( signals.length , 0.0 ); 
    for( int i = 0; i < signals.size(); i++ )
    {
        signal_values[i] = SignalBehavior.getSingleSignal( pCell, signals.get(i) );
    }

        // now, get live/dead value 
        boolean dead = SignalBehavior.getSingleSignal( pCell, "dead" ) != 0;

        double out = evaluate( signal_values , dead ); 

        // new March 27, 2023 
        // if the rule was found to not apply, then just get the prior value 
        if( out < -9e90 )
        {
            out = SignalBehavior.getSinglBehavior( pCell, this.behavior );
        }

        return out; 
    }

    void apply(Cell pCell) throws Exception
    {
        // evaluate the rule 
        double param = evaluate( pCell );

        // apply it ot the appropriate behavior 
        SignalBehavior.setSingleBehavior( pCell, behavior, param );
    }

    void sync_to_CellDefinition(CellDefinition pCD)
    {
        if( pCD == null )
            return;
        cell_type = pCD.name;
        //        pCellDefinition = pCD;

        // sync base behavior 
        base_value = SignalBehavior.get_single_base_behavior( pCD, behavior );
    }

    void sync_to_CellDefinition(String cell_name)
    {
        sync_to_CellDefinition( CellDefinition.getCellDefinition( cell_name ) );
    }

    int find_signal(String name)
    {
        if( signals_map.containsKey( name ) )
            signals_map.get( name );
        return -1;
        //        auto search = signals_map.find( name );
        //
        //        if( search == signals_map.end() )
        //        {
        //            return -1;
        //        }
        //
        //        return search.second;
    }

    void set_half_max(String name, double hm)
    {
        int n = find_signal( name );
        if( n < 0 )
        {
            return;
        }

        half_maxes.set( n, hm );//[n] = hm;

        if( responses.get( n ) == true )
        {
            for( int m = 0; m < up_signals.size(); m++ )
            {
                if( up_signals.get( m ).equals( name ) )
                {
                    up_half_maxes.set( m, hm );//[m] = hm;
                }
            }
        }
        else
        {
            for( int m = 0; m < down_signals.size(); m++ )
            {
                if( down_signals.get( m ).equals( name ) )
                {
                    down_half_maxes.set( m, hm );//[m] = hm;
                }
            }
        }
        return;
    }

    void set_hill_power(String name, double hp)
    {
        int n = find_signal( name );
        if( n < 0 )
        {
            return;
        }

        hill_powers.set( n, hp );//[n] = hp;
        if( responses.get( n ) == true )
        {
            for( int m = 0; m < up_signals.size(); m++ )
            {
                if( up_signals.get( m ).equals( name ) )//[m] == name )
                {
                    up_hill_powers.set( m, hp );//[m] = hp;
                }
            }
        }
        else
        {
            for( int m = 0; m < down_signals.size(); m++ )
            {
                if( down_signals.get( m ).equals( name ) )//[m] == name )
                {
                    down_hill_powers.set( m, hp );//[m] = hp;
                }
            }
        }
        return;
    }

    void set_response(String name, String response)
    {
        int n = find_signal( name );
        if( n < 0 )
        {
            return;
        }

        boolean bResponse = false; // true if up-regulate, false if down
        if( response == "increase" || response == "increases" || response == "promotes" )
        {
            bResponse = true;
        }

        // this is already my response? if so exit 
        if( bResponse == responses.get( n ) )
        {
            return;
        }

        if( responses.get( n ) )//[n] == true )
        {
            // need to switch from up to down 
            // find current index 
            int ci = -1;
            for( int m = 0; m < up_signals.size(); m++ )
            {
                if( up_signals.get( m ).equals( name ) )//[m] == name )
                {
                    ci = m;
                }
            }
            // swap last inot that position 
            up_half_maxes.set( ci, up_half_maxes.get( up_half_maxes.size() - 1 ) );//.back() );
            up_hill_powers.set( ci, up_hill_powers.get( up_hill_powers.size() - 1 ) );
            up_signals.set( ci, up_signals.get( up_signals.size() - 1 ) );
            up_applies_to_dead_cells.set( ci, up_applies_to_dead_cells.get( up_applies_to_dead_cells.size() - 1 ) );

            //                        up_half_maxes[ci] = up_half_maxes.back();
            //            up_hill_powers[ci] = up_hill_powers.back();
            //            up_signals[ci] = up_signals.back();
            //            up_applies_to_dead_cells[ci] = up_applies_to_dead_cells.back();

            // reduce size by one

            up_half_maxes.remove( up_half_maxes.size() - 1 );//.pop_back();
            up_hill_powers.remove( up_hill_powers.size() - 1 );//pop_back();
            up_signals.remove( up_signals.size() - 1 );//pop_back();
            up_applies_to_dead_cells.remove( up_applies_to_dead_cells.size() - 1 );//pop_back();

            // move to the other side 

            down_half_maxes.add( half_maxes.get( n ) );//[n] );
            down_hill_powers.add( hill_powers.get( n ) );
            down_signals.add( signals.get( n ) );
            down_applies_to_dead_cells.add( applies_to_dead_cells.get( n ) );
        }
        else
        {
            // need to switch from down to up 
            // find current index 
            int ci = -1;
            for( int m = 0; m < down_signals.size(); m++ )
            {
                if( down_signals.get( m ).equals( name ) )//[m] == name )
                {
                    ci = m;
                }
            }
            // swap last inot that position 
            down_half_maxes.set( ci, down_half_maxes.get( down_half_maxes.size() - 1 ) );//.back() );
            down_hill_powers.set( ci, down_hill_powers.get( down_hill_powers.size() - 1 ) );
            down_signals.set( ci, down_signals.get( down_signals.size() - 1 ) );
            down_applies_to_dead_cells.set( ci, down_applies_to_dead_cells.get( down_applies_to_dead_cells.size() - 1 ) );
            //            down_half_maxes[ci] = down_half_maxes.back();
            //            down_hill_powers[ci] = down_hill_powers.back();
            //            down_signals[ci] = down_signals.back();
            //            down_applies_to_dead_cells[ci] = up_applies_to_dead_cells.back();

            // reduce size by one
            down_half_maxes.remove( down_half_maxes.size() - 1 );//.pop_back();
            down_hill_powers.remove( down_hill_powers.size() - 1 );//pop_back();
            down_signals.remove( down_signals.size() - 1 );//pop_back();
            down_applies_to_dead_cells.remove( down_applies_to_dead_cells.size() - 1 );//pop_back();
            //            down_half_maxes.pop_back();
            //            down_hill_powers.pop_back();
            //            down_signals.pop_back();
            //            down_applies_to_dead_cells.pop_back();

            // move to the other side 

            up_half_maxes.add( half_maxes.get( n ) );
            up_hill_powers.add( hill_powers.get( n ) );
            up_signals.add( signals.get( n ) );
            up_applies_to_dead_cells.add( applies_to_dead_cells.get( n ) );
        }

        responses.set( n, bResponse );// [n] = bResponse;
    }
}
