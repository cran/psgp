/***************************************************************************
 *   AstonGeostats, algorithms for low-rank geostatistical models          *
 *                                                                         *
 *   Copyright (C) Ben Ingram, 2008-2009                                   *
 *                                                                         *
 *   Ben Ingram, IngramBR@Aston.ac.uk                                      *
 *   Neural Computing Research Group,                                      *
 *   Aston University,                                                     *
 *   Aston Street, Aston Triangle,                                         *
 *   Birmingham. B4 7ET.                                                   *
 *   United Kingdom                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef EXPRESSIONNODE_H_
#define EXPRESSIONNODE_H_

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

enum Operator {NOP, PLUS, MINUS, TIMES, DIVIDE, POWER, UNARYMINUS, SIN, COS, EXP, LOG};
enum Nodetype {LITERAL, VARIABLE, OPERATOR};


using namespace std;

class ExpressionNode
{
public:
	
	ExpressionNode(string id);
	ExpressionNode(double value);
	ExpressionNode(Operator opt, ExpressionNode* lval, ExpressionNode* rval);
	virtual ~ExpressionNode();

	double evaluate(double var);
	string toString();

private:

	double evaluateNode();
	double evaluateOperatorNode();
	double evaluateValueNode();
	double evaluateVariableNode();

	string stringNode();
	string stringOperatorNode();
	string stringValueNode();
	string stringVariableNode();

	bool isOperator()
	{
		return nodeType == OPERATOR;
	}

	bool isLeaf()
	{
		return nodeType != OPERATOR;
	}

	bool isLiteral()
	{
		return nodeType == LITERAL;
	}
      
	bool isVariable()
	{
		return nodeType == VARIABLE;
	}
    
	Nodetype getNodeType()
	{
		return nodeType;
	}
    
	Operator getOperator()
	{
		return operationType;
	}
    
	double getDoubleValue()
	{
		return literalValue;
	}
    
	string getVarName()
	{
		return variableName;
	}

	ExpressionNode* leftLeaf;
 	ExpressionNode* rightLeaf;
	Nodetype nodeType;
	Operator operationType;
    
	double literalValue;
	string variableName;

	static double variable;
};

#endif

