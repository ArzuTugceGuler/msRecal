/* 
 * File:   pepXMLDelegates.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 14:45
 */

#ifndef PEPXMLDELEGATES_H
#define	PEPXMLDELEGATES_H

typedef struct pepxml_delegate
{
	char* key;
	char* context;
	void* (*pread_function)(char*);
	char* (*write_function)(void*);
	void* (*property_function)(char*, void*);

} delegate_type, *pdelegate_type;


/* List of delegates */
typedef struct pepxml_delegate_list
{
	pdelegate_type*	dlist;	/* list of pointers to delegates */
	int				size;	/* size of the list */

} delegate_list, *pdelegate_list;


/* Function that creates and configures a delegate */
pdelegate_type make_delegate(char* key, char* context, void* (*pread_function)(char*), char* (*write_function)(void*), void* (*property_function)(char*, void*));


/* Function that creates a delegate list */
pdelegate_list make_delegate_list();

/* Function that adds a delegate to the list */
void add_delegate_to_list(pdelegate_list list, pdelegate_type delegate_t);

/* Function that helps finding the correct delegate, NULL if key is not found in the array */
pdelegate_type get_delegate(char* key, char* context, pdelegate_list list);

#endif	/* PEPXMLDELEGATES_H */

