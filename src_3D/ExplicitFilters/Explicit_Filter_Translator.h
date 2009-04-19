/*
 *  Explicit_Filter_Translator.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _EXPLICIT_FILTER_TRANSLATOR_INCLUDED
#define _EXPLICIT_FILTER_TRANSLATOR_INCLUDED


#define TRANSLATE_MEMBER 0
#define TRANSLATE_ALL    1

template<class ObjectType, class Member_Pointer>
class Explicit_Filter_Translator {
private:
    int N;
    int type;
    ObjectType Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		/*!< pointer to the class member function */
    
public:
    Explicit_Filter_Translator(void);
    Explicit_Filter_Translator(ObjectType object, Member_Pointer member){
        Set(object,member);
    }
    
    Explicit_Filter_Translator(ObjectType object) { 
        Set(object) ;
    }
    void Set(ObjectType object, Member_Pointer member) {
        Obj = object;
        Ptr = member;
        N = 1;
        type = TRANSLATE_MEMBER;
    }
    void Set(ObjectType object) {
        N = object.NumVar() ; type = TRANSLATE_ALL;
    }
    
    RowVector operator() (Cell3D &theCell) {
        return Obj[theCell.I][theCell.J][theCell.K].*Ptr;
    }
    
    RowVector operator() ( int i,  int j,  int k) {
        if (type == TRANSLATE_MEMBER) {
            RowVector temp(N);
            temp(0) = (Obj[i][j][k]).*Ptr;
            return temp;
        } else {
            RowVector temp(N);
            /*for (int n=0; n<N; n++) {
                temp(n) = Obj[i][j][k][n];
            }*/
            cout << "not supported" << endl;
            return temp;
        }
    }
};


#endif
