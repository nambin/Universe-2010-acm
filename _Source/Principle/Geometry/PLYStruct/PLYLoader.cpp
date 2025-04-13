#include <stdafx.h>
#include <Principle/Geometry/PLYLibrary/PLYLoader.h>

namespace plyloader {

void plyLoader(FILE *fp, int &nverts, int &nfaces, Vertex **&vlist, Face **&flist, int &per_vertex_color, int &has_normals) {

	// Read in the original PLY object 
	PlyFile *in_ply = read_ply(fp);

	int elem_count;
	for (int i = 0; i < in_ply->num_elem_types; ++i) {

		// prepare to read the i'th list of elements 
		char *elem_name = setup_element_read_ply (in_ply, i, &elem_count);

		if (equal_strings ("vertex", elem_name)) {

			// create a vertex list to hold all the vertices 
			vlist = (Vertex **) malloc (sizeof (Vertex *) * elem_count);
			nverts = elem_count;

			// set up for getting vertex elements 
			setup_property_ply (in_ply, &vert_props[0]);
			setup_property_ply (in_ply, &vert_props[1]);
			setup_property_ply (in_ply, &vert_props[2]);

			for (int j = 0; j < in_ply->elems[i]->nprops; ++j) {
				PlyProperty *prop = in_ply->elems[i]->props[j];				
				if (equal_strings ("red", prop->name) || equal_strings ("r", prop->name)) {
					setup_property_ply (in_ply, &vert_props[3]);
					per_vertex_color = 1;
				}
				else if (equal_strings ("green", prop->name) || equal_strings ("g", prop->name)) {
					setup_property_ply (in_ply, &vert_props[4]);
					per_vertex_color = 1;
				}
				else if (equal_strings ("blue", prop->name) || equal_strings ("b", prop->name)) {
					setup_property_ply (in_ply, &vert_props[5]);
					per_vertex_color = 1;
				}
				else if (equal_strings ("nx", prop->name)) {
					setup_property_ply (in_ply, &vert_props[6]);
					has_normals = 1;
				}
				else if (equal_strings ("ny", prop->name)) {
					setup_property_ply (in_ply, &vert_props[7]);
					has_normals = 1;
				}
				else if (equal_strings ("nz", prop->name)) {
					setup_property_ply (in_ply, &vert_props[8]);
					has_normals = 1;
				}
				else if (equal_strings ("confidence", prop->name)) {
					setup_property_ply (in_ply, &vert_props[9]);
					per_vertex_color = 1;
				}
			}
			//vert_other = get_other_properties_ply (in_ply, offsetof(Vertex, other_props));

			// grab all the vertex elements 
			for (int j = 0; j < elem_count; j++) {
				vlist[j] = (Vertex *) malloc (sizeof (Vertex));
				get_element_ply (in_ply, (void *) vlist[j]);
			}			
		}
		else if (equal_strings ("face", elem_name)) {

			// create a list to hold all the face elements 
			flist = (Face **) malloc (sizeof (Face *) * elem_count);
			nfaces = elem_count;

			// set up for getting face elements 
			setup_property_ply (in_ply, &face_props[0]);

			//face_other = get_other_properties_ply (in_ply, offsetof(Face, other_props));

			// grab all the face elements 
			for (int j = 0; j < elem_count; j++) {
				flist[j] = (Face *) malloc (sizeof (Face));
				get_element_ply (in_ply, (void *) flist[j]);
			}
		}
		//else
			//get_other_element_ply (in_ply);
	}

	// Swap bits if needed
	if(in_ply->file_type == PLY_BINARY_BE) {
		for(int i=0;i<nverts;++i) {
			basis::BYTE_SWAP(vlist[i]->x);basis::BYTE_SWAP(vlist[i]->y);basis::BYTE_SWAP(vlist[i]->z);
			basis::BYTE_SWAP(vlist[i]->confidence);
			basis::BYTE_SWAP(vlist[i]->nx);basis::BYTE_SWAP(vlist[i]->ny);basis::BYTE_SWAP(vlist[i]->nz);
		}
		for(int i=0;i<nfaces;++i) {
			for(int j=0;j<flist[i]->nverts;++j) {
				basis::BYTE_SWAP(flist[i]->verts[j]);
			}
		}
	}

	close_ply (in_ply);
	free_ply (in_ply);
}

std::ostream &operator<<(std::ostream &os, const Vertex &v) {
	os << "XYZ : " << v.x << " " << v.y << " " << v.z << std::endl;
	os << "RGB : " << (int)v.r << " " << (int)v.g << " " << (int)v.b << std::endl;
	os << "NORMAL : " << v.nx << " " << v.ny << " " << v.nz << std::endl;
	os << "CONF : " << v.confidence << std::endl;

	return os;
}

std::ostream &operator<<(std::ostream &os, const Face &v) {
	os << "nVerts : " << v.nverts << " _ " << "vertIdx : ";
	for(int i=0;i<v.nverts;++i) 
		os << v.verts[i] << " ";
	os << std::endl;
	return os;
}

};		// plyloader
