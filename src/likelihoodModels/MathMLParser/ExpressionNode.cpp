#include "ExpressionNode.h"

#define BAD_NUMBER 1.0
	
ExpressionNode::ExpressionNode(string id) : leftLeaf(NULL), rightLeaf(NULL), nodeType(VARIABLE), variableName(id)
{

}

ExpressionNode::ExpressionNode(double value) : leftLeaf(NULL), rightLeaf(NULL), nodeType(LITERAL), literalValue(value)
{

}

ExpressionNode::ExpressionNode(Operator opt, ExpressionNode* lval, ExpressionNode* rval) : leftLeaf(lval), rightLeaf(rval), nodeType(OPERATOR) , operationType(opt)
{

}

ExpressionNode::~ExpressionNode()
{
	if(leftLeaf != NULL)
	{
		delete leftLeaf;
	}
	if(rightLeaf != NULL)
	{
		delete rightLeaf;
	}
}



double ExpressionNode::evaluate(double var)
{
	ExpressionNode::variable = var;
	double result = evaluateNode();
	return result;
}

double ExpressionNode::evaluateNode()
{
	switch (nodeType)
	{
		case LITERAL : 
			return(this->evaluateValueNode());

		case VARIABLE : 
			return(this->evaluateVariableNode());

		case OPERATOR :
			return(this->evaluateOperatorNode());

		default :
			cerr << "ExpressionNode: Node type not implemented" << endl;
	}
	return BAD_NUMBER;
}

double ExpressionNode::evaluateOperatorNode()
{
	double res;
       
	switch (operationType)
	{
		case PLUS:
			res = leftLeaf->evaluateNode() + rightLeaf->evaluateNode();
			break;

		case MINUS:
			res = leftLeaf->evaluateNode() - rightLeaf->evaluateNode();
			break;
              
		case TIMES:
			res = leftLeaf->evaluateNode() * rightLeaf->evaluateNode();
			break;

		case DIVIDE:
			res = leftLeaf->evaluateNode() / rightLeaf->evaluateNode();			
			break;

		case POWER:
			res = pow(leftLeaf->evaluateNode(), rightLeaf->evaluateNode());			
			break;

		case EXP:
			res = exp(leftLeaf->evaluateNode());
			break;

		case LOG:
			res = log(leftLeaf->evaluateNode());
			break;

		case SIN:
			res = sin(leftLeaf->evaluateNode());
			break;

		case COS:
			res = cos(leftLeaf->evaluateNode());
			break;

		case UNARYMINUS:
			res = -(leftLeaf->evaluateNode());
			break;

		case NOP:
			res = BAD_NUMBER;
			cerr << "ExpressionNode: NOP" << endl;
			break;

		default:
			res = BAD_NUMBER;
			cerr << "ExpressionNode: Unknown operator" << endl;
	}

	return res;
}

string ExpressionNode::stringNode()
{
	switch (nodeType)
	{
		case LITERAL : 
			return(this->stringValueNode());

		case VARIABLE : 
			return(this->stringVariableNode());

		case OPERATOR :
			return(this->stringOperatorNode());

		default :
			cerr << "ExpressionNode: Node type not implemented" << endl;
	}
	return "***ERROR***";
}


string ExpressionNode::toString()
{
	return stringNode();
}


string ExpressionNode::stringOperatorNode()
{
	string str;

	switch (operationType)
	{
		case PLUS:
			str = "(" + leftLeaf->stringNode() + "+" + rightLeaf->stringNode() + ")";
			break;

		case MINUS:
			str = "(" + leftLeaf->stringNode() + "-" + rightLeaf->stringNode() + ")";
			break;
              
		case TIMES:
			str = "(" + leftLeaf->stringNode() + "*" + rightLeaf->stringNode() + ")";
			break;

		case DIVIDE:
			str = "(" + leftLeaf->stringNode() + "/" + rightLeaf->stringNode() + ")";			
			break;

		case POWER:
			str = "pow(" + leftLeaf->stringNode() + ", " + rightLeaf->stringNode() + ")";			
			break;

		case EXP:
			str = "exp(" + leftLeaf->stringNode() + ")";
			break;

		case LOG:
			str = "log(" + leftLeaf->stringNode() + ")";
			break;

		case SIN:
			str = "sin(" + leftLeaf->stringNode() + ")";
			break;

		case COS:
			str = "cos(" + leftLeaf->stringNode() + ")";
			break;

		case UNARYMINUS:
			str = "(-" + leftLeaf->stringNode() + ")";
			break;

		case NOP:
			str = BAD_NUMBER;
			cerr << "ExpressionNode: NOP" << endl;
			break;

		default:
			str = BAD_NUMBER;
			cerr << "ExpressionNode: Unknown operator" << endl;
	}

	return str;
}

string ExpressionNode::stringValueNode()
{
	std::string s;
	std::ostringstream ss;
	ss << literalValue;
	s = ss.str();
	return s;
}

string ExpressionNode::stringVariableNode()
{
	return variableName;
}

double ExpressionNode::evaluateValueNode()
{
	return literalValue;
}

double ExpressionNode::evaluateVariableNode()
{
	return ExpressionNode::variable;
}

double ExpressionNode::variable = 0.0;
