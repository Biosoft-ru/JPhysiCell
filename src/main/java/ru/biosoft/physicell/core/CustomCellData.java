package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.HashMap;
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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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
public class CustomCellData
{
    //    std::unordered_map<std::string,int> name_to_index_map; 
    Map<String, Integer> name_to_index_map = new HashMap<>();

    /// /  std::unordered_map<std::string,int> vector_name_to_index_map; 

    //    friend std::ostream&operator<<(std::ostream&os,const Custom_Cell_Data&ccd); // done 
    //    public:std::vector<Variable> variables;
    //    std::vector<Vector_Variable> vector_variables;

    //    Variable[] variables;
    //    VectorVariable[] vector_variables;
    List<Variable> variables = new ArrayList<>();
    List<VectorVariable> vector_variables = new ArrayList<VectorVariable>();

    public CustomCellData()
    {
    //  std::cout << __FUNCTION__ << "(default)" << std::endl; 
    //        variables = new Variable[0];//.resize(0); 
    //        vector_variables = new VectorVariable[0];//.resize(0); 
        
        name_to_index_map.clear(); 
    //  vector_name_to_index_map.clear();
    }

    public CustomCellData(CustomCellData ccd)
    {
    //  std::cout << __FUNCTION__ << "(copy)" << std::endl; 
        variables = ccd.variables; 
        vector_variables = ccd.vector_variables; 
        name_to_index_map= ccd.name_to_index_map; 
    }

    int add_variable(Variable v)
    {
        int n = variables.size(); 
        variables.add( v );//.push_back( v ); 
        name_to_index_map.put( v.name, n );//[ v.name ] = n; 
        return n; 
    }

    int add_variable(String name, String units, double value)
    {
        int n = variables.size(); 
        Variable variable = new Variable();
        variable.name = name;
        variable.units = units;
        variable.value = value;
        variables.add( variable );
        name_to_index_map.put( name, n );
        //        variables.resize( n+1 ); 
        //        variables[n].name = name; 
        //        variables[n].units = units; 
        //        variables[n].value = value; 
        //        name_to_index_map[ name ] = n; 
        return n; 
    }

    int add_variable(String name , double value )
    {
        //        int 
        return add_variable( name, "dimensionless", value );
        //        variables.resize( n+1 ); 
        //        variables[n].name = name; 
        //        variables[n].units = "dimensionless"; 
        //        variables[n].value = value; 
        //        name_to_index_map[ name ] = n; 
        //        return n; 
    }

    int add_vector_variable( VectorVariable v )
    {
        int n = vector_variables.size();
        vector_variables.add( v );
        name_to_index_map.put( v.name, n );
        return n;
        //        int n = vector_variables.size(); 
        //        vector_variables.push_back( v ); 
    //  vector_name_to_index_map[ v.name ] = n; 
    //        return n; 
    }

    int add_vector_variable( String name , String units , double[] value )
    {
        VectorVariable v = new VectorVariable();
        v.name = name;
        v.units = units;
        v.value = value.clone();
        return add_vector_variable( v );
        //        int n = vector_variables.size(); 
        //        vector_variables.resize( n+1 ); 
        //        vector_variables[n].name = name; 
        //        vector_variables[n].units = units; 
        //        vector_variables[n].value = value; 
    //  vector_name_to_index_map[ name ] = n; 
        //        return n; 
    }

    int add_vector_variable( String name , double[] value )
    {
        return add_vector_variable( name, "dimensionless", value );
        //        int n = vector_variables.size(); 
        //        vector_variables.resize( n+1 ); 
        //        vector_variables[n].name = name; 
        //        vector_variables[n].units = "dimensionless"; 
        //        vector_variables[n].value = value; 
        //    //  vector_name_to_index_map[ name ] = n; 
        //        return n; 
    }

    int find_variable_index( String name )
    {
        if( name_to_index_map.containsKey( name ) )
            return name_to_index_map.get( name );
        return -1;
        // this should return -1 if not found, not zero 
        //        auto out = name_to_index_map.find( name ); 
        //        if( out != name_to_index_map.end() )
        //        { return out->second; }
        //        return -1; 
    }

    /*
    int Custom_Cell_Data::find_vector_variable_index( std::string name )
    {
        return vector_name_to_index_map[ name ]; 
    }
    */

    int find_vector_variable_index( String name )
    {
        for( int i = 0; i < vector_variables.size(); i++ )
        {
            if( vector_variables.get( i ).name.equals( name ) )
                return i;
        }
        return -1;
        //        int n = 0; 
        //        while( n < vector_variables.size() )
        //        {
        //            if( std::strcmp( vector_variables[n].name.c_str() , name.c_str() ) == 0 )
        //            { return n; } 
        //            n++; 
        //        }
        //        
        //        return -1; 
    }

    public double get(int i)
    {
        return variables.get( i ).value;
    }
    public void set(int i, double val)
    {
        variables.get( i ).value = val;
    }
    
    public void set(String name, double val)
    {
        int index = name_to_index_map.get( name );
        variables.get( index ).value = val;
    }
    
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Custom data (scalar): "); 
        for( int i = 0; i < variables.size(); i++ )
        {
            sb.append( i + ": " + variables.get( i ) );
        }

        sb.append( "Custom data (vector): "); 
        for( int i = 0; i < vector_variables.size(); i++ )
        {
            sb.append( i + ": " + vector_variables.get( i ) );
        }
        
        return sb.toString();
    }

    public static class Variable
    {
        String name;
        double value;
        String units;
        boolean conserved_quantity;

        public Variable()
        {
            name = "unnamed";
            units = "dimensionless";
            value = 0.0;
            conserved_quantity = false;
            return;
        }

        @Override
        public String toString()
        {
            return name + ": " + value + " " + units; 
        }
    }

    public class VectorVariable
    {
        String name;
        double[] value;
        String units;
        boolean conserved_quantity;

        public VectorVariable()
        {
            name = "unnamed";
            units = "dimensionless";
            value = new double[3];//.resize(3, 0.0 );
            conserved_quantity = false;
            return;
        }

        public String toString()
        {
            StringBuilder sb = new StringBuilder();
            sb.append( name + ": ");
            if( value.length == 0 )
            { sb.append( "[empty]");
            return sb.toString();
            }
            
            /*
            if( v.value.size() == 1 )
            { os << v.value[0] << " (" << v.units << ")"; return os;  }
        */
            for( int i=0; i < value.length-1 ; i++ )
            { 
                sb.append( value[i]+"," );
//                os << v.value[i] << ","; 
            }
            sb.append( value[value.length - 1] + " (" + units + " )" );
//            os << v.value[v.value.size()-1] << " (" << v.units << ")"; 
            return sb.toString(); 
        }
    }
}