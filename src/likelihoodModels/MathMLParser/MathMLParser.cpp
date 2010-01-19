
#include "MathMLParser.h"

MathMLParser::MathMLParser()
{

}

MathMLParser::~MathMLParser()
{

}

bool MathMLParser::isUnaryOperator(char *op)
{
	bool found = false;
	int entries = 6;
	string Unary_Operators[6] = {"abs", "exp", "log", "minus", "cos", "sin"};
	string name(op);

	for(int i = 0 ; i < entries ; i++)
	{
		if(!Unary_Operators[i].compare(name))
		{
			found = true;
		}
	}
	return(found);
}


bool MathMLParser::isBinaryOperator(char *op)
{
	bool found = false;
	int entries = 5;
	string Binary_Operators[5] = {"plus", "minus", "times", "divide", "power"}; 
	string name(op);

	for(int i = 0 ; i < entries ; i++)
	{
		if(!Binary_Operators[i].compare(name))
		{
			found = true;
		}
	}
	return(found);
}

bool MathMLParser::isOperator(char *name)
{
	if (isBinaryOperator(name))
	{
		return true;
	}
	if (isUnaryOperator(name))
	{
		return true;
	}
	return false;
}


DOMNode* MathMLParser::skipJunkNodes(DOMNode *node)
{
	char *name;
	bool found = false;
	int count;
	DOMNode *n = node;
	DOMText *t;
	while ((n != NULL) && (found == false))
	{
		switch (n->getNodeType())
		{
			case DOMNode::ATTRIBUTE_NODE:
				break;
			case DOMNode::CDATA_SECTION_NODE:
				break;
			case DOMNode::COMMENT_NODE:
				break;
			case DOMNode::DOCUMENT_FRAGMENT_NODE:
				break;
			case DOMNode::DOCUMENT_TYPE_NODE:
				break;
			case DOMNode::ELEMENT_NODE:
				found = true;
				break;
			case DOMNode::NOTATION_NODE:
				break;
			case DOMNode::TEXT_NODE:
				t = (DOMText *)n;
				name = XMLString::transcode(t->getData());
				count = 0;
				for(int i = 0 ; i < (int)strlen(name) ; i++)
				{
					if((name[i] == ' ') || (name[i] == '\n'))
					{
						count++;
					}
				}
				if((strlen(name) - count) > 0)
				{
					found = true;
				}
				break;
			case DOMNode::PROCESSING_INSTRUCTION_NODE:
				break;
			default:
				break;
		}

		if (found == false)
		{
			n = n->getNextSibling();
		}

	}
	return n;
}


DOMNode* MathMLParser::skipTextNodes(DOMNode *node)
{
	DOMNode *n = node;
	if(n == NULL)
	{
		return NULL;
	}

	while(true)
	{
		n = skipJunkNodes(n);
		if((n != NULL) && (n->getNodeType() == DOMNode::TEXT_NODE))
		{
			n = n->getNextSibling();
		}
		else
		{
			break;
		}
	}
	return n;
}

DOMNode* MathMLParser::nextElementNode(DOMNode *node)
{
	if (node == NULL)
	{
		return NULL;
	}
	return skipTextNodes(node->getNextSibling());
}

ExpressionNode* MathMLParser::parseIdentifier(DOMNode *node)
{

	DOMNode *child = node->getFirstChild();
	child = skipJunkNodes(child);
        
	if (child == NULL)
	{
		cerr << "MathMLParser::parseIdentifier: empty identifier element <ci>" << endl;
	}

	if (child->getNodeType() != DOMNode::TEXT_NODE)
	{
		cerr << "MathMLParser::parseIdentifier: text node (identifier) expected" << endl;
	}

	string str(XMLString::transcode(((DOMText *)child)->getData()));
	if (str.length() == 0)
	{
		cerr << "MathMLParser::parseIdentifier: whitespace is not an allowed identifier (ci)" << endl;
	}
        
	return new ExpressionNode(str);
}

ExpressionNode* MathMLParser::parseNumber(DOMNode *node)
{
	DOMNode *child = NULL;
	double dblval;

	NumberType type = REAL;

	ExpressionNode *value = NULL;
        
	child = node->getFirstChild();
	child = skipJunkNodes(child);
        
	if (child == NULL)
	{
		cerr << "MathML::parsenumber: Child==NULL" << endl;
	}

	DOMNamedNodeMap *nnm = node->getAttributes();
	DOMAttr *typeAttr = (DOMAttr *)nnm->getNamedItem(XMLString::transcode("type"));
	DOMAttr *baseAttr = (DOMAttr *)nnm->getNamedItem(XMLString::transcode("base"));

	if (baseAttr != NULL)
	{
		type = INTEGER;
	}
	else
	{
		if(typeAttr == NULL)
			type = REAL;
		else if(!strcmp(XMLString::transcode(typeAttr->getValue()),"integer"))
			type = INTEGER;
		else if(!strcmp(XMLString::transcode(typeAttr->getValue()),"real"))
			type = REAL;
		else if(!strcmp(XMLString::transcode(typeAttr->getValue()),"rational"))
			type = RATIONAL;
		else if(!strcmp(XMLString::transcode(typeAttr->getValue()),"e-notation"))
			type = SCIENTIFIC;
      else
			type = UNKNOWN;
	}

	switch (type)
	{
		case REAL:
			dblval = strtod(XMLString::transcode(((DOMText *)child)->getData()), NULL);
			value = new ExpressionNode(dblval);
			break;
		case INTEGER:
		case RATIONAL:
		case SCIENTIFIC:
		default:
			cerr << "MathML::parseNumber: Unsupported number type: " << type << endl;
	}
	return value;
}


Operator MathMLParser::parseOperator(DOMNode *node)
{
	Operator op;
	char *operationCode = XMLString::transcode(node->getNodeName());

	// wouldn't switch on strings be handy here?
	if(!strcmp(operationCode, "plus"))
		op = PLUS;
	else if(!strcmp(operationCode, "minus"))
		op = MINUS;
	else if(!strcmp(operationCode, "times"))
		op = TIMES;
	else if(!strcmp(operationCode, "divide"))
		op = DIVIDE;
	else if(!strcmp(operationCode, "power"))
		op = POWER;
	else if(!strcmp(operationCode, "unaryminus"))
		op = UNARYMINUS;
	else if(!strcmp(operationCode, "sin"))
		op = SIN;
	else if(!strcmp(operationCode, "cos"))
		op = COS;
	else if(!strcmp(operationCode, "exp"))
		op = EXP;
	else if(!strcmp(operationCode, "log"))
		op = LOG;
	else
	{
		op = NOP;
		cerr << "MathML::parseOperator: Unknown operator" << endl;
	}
	return op;
} 

ExpressionNode* MathMLParser::parseApply(DOMNode *node)
{
	char *opName;
	DOMNode *child, *opNode;
	ExpressionNode *leftChild, *rightChild, *rightExpression;
	Operator op;

	child = node->getFirstChild();
	child = skipJunkNodes(child);       

	if (child == NULL)
	{
		cerr << "MathML::parseApply: operator node expected - Chile==NULL" << endl;
	}

	if(child->getNodeType() != DOMNode::ELEMENT_NODE)
	{
		cerr << "MathML::parseApply: operator node expected - Not ELEMENT_NODE" << endl;
	}

	if(isOperator(XMLString::transcode(child->getNodeName())) == false)
	{
		cerr << "MathML::parseApply: operator node expected - Not Operator" << endl;
	}

	op = parseOperator(child);

	opNode = child;
	opName = XMLString::transcode(opNode->getNodeName());
	child = nextElementNode(child);
        
	if((child == NULL) || (child->getNodeType() != DOMNode::ELEMENT_NODE))
	{
		cerr << "MathML::parseApply: Left-value expected" << endl;
	}

	leftChild = parseExpression(child);
	rightChild = NULL;
	child = nextElementNode(child);

	while(child != NULL)
	{
		if(child->getNodeType() != DOMNode::ELEMENT_NODE)
		{
			cerr << "MathML::parseApply: Wrong element type" << endl;
		}
            
		rightExpression = parseExpression(child);
            
		if(rightChild == NULL)
		{
			if(isUnaryOperator(opName) && !isBinaryOperator(opName))
			{
				cerr << "MathML::parseApply: Too many arguments to unary operator." << endl;
			}
			rightChild = rightExpression;
			rightExpression = NULL;
		}
		else
		{
			if (isBinaryOperator(opName))
			{
				cerr << "MathML::parseApply: Too many arguments to binary operator." << endl;
			}
			rightChild = new ExpressionNode(op, rightChild, rightExpression);
			rightExpression = NULL;
		}
		child = nextElementNode(child);
	}

	return new ExpressionNode(op, leftChild, rightChild);
}

ExpressionNode* MathMLParser::parseExpression(DOMNode *node)
{
	char *name = XMLString::transcode(node->getNodeName());

	ExpressionNode *expr;

	// need to perform case-insensitive string comparison        
	if(!strcmp(name, "ci"))
	{
		expr = parseIdentifier(node);
	}
	else
	{
		if(!strcmp(name, "cn"))
		{
			expr = parseNumber(node);
		}
		else
		{
			if(!strcmp(name, "apply"))
			{
				expr = parseApply(node);
			}
			else
			{
				expr = NULL;
				cerr << "MathMLParser::parseExpression: ci, cn or apply expected" << endl;
			}
		}
	}
	return expr;
}

