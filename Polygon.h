// Polygon.h: interface for the Polygon class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYGON_H__6B661CC4_2E36_4F4F_8F4D_872D408901D6__INCLUDED_)
#define AFX_POLYGON_H__6B661CC4_2E36_4F4F_8F4D_872D408901D6__INCLUDED_



#define GPC_EPSILON (DBL_EPSILON)

#define GPC_VERSION "2.32"

;
//##ModelId=4E4A89B301B9
typedef enum                        /* Set operation type                */
{
  GPC_DIFF,                         /* Difference                        */
  GPC_INT,                          /* Intersection                      */
  GPC_XOR,                          /* Exclusive or                      */
  GPC_UNION                         /* Union                             */
} gpc_op;

//##ModelId=4E4A89B301C5
typedef struct                      /* Polygon vertex structure          */
{
  double              x;            /* Vertex x component                */
  double              y;            /* vertex y component                */
} gpc_vertex;

//##ModelId=4E4A89B301C7
typedef struct                      /* Vertex list structure             */
{
  int                 num_vertices; /* Number of vertices in list        */
  gpc_vertex         *vertex;       /* Vertex array pointer              */
} gpc_vertex_list;

//##ModelId=4E4A89B301D4
typedef struct                      /* Polygon set structure             */
{
  int                 num_contours; /* Number of contours in polygon     */
  int                *hole;         /* Hole / external contour flags     */
  gpc_vertex_list    *contour;      /* Contour array pointer             */
} gpc_polygon;

//##ModelId=4E4A89B301D6
typedef struct                      /* Tristrip set structure            */
{
  int                 num_strips;   /* Number of tristrips               */
  gpc_vertex_list    *strip;        /* Tristrip array pointer            */
} gpc_tristrip;


/*
===========================================================================
                       Public Function Prototypes
===========================================================================
*/

void gpc_read_polygon        (FILE            *infile_ptr, 
                              int              read_hole_flags,
                              gpc_polygon     *polygon);

void gpc_write_polygon       (FILE            *outfile_ptr,
                              int              write_hole_flags,
                              gpc_polygon     *polygon);

void gpc_add_contour         (gpc_polygon     *polygon,
                              gpc_vertex_list *contour,
                              int              hole);

void gpc_polygon_clip        (gpc_op           set_operation,
                              gpc_polygon     *subject_polygon,
                              gpc_polygon     *clip_polygon,
                              gpc_polygon     *result_polygon);

void gpc_tristrip_clip       (gpc_op           set_operation,
                              gpc_polygon     *subject_polygon,
                              gpc_polygon     *clip_polygon,
                              gpc_tristrip    *result_tristrip);

void gpc_polygon_to_tristrip (gpc_polygon     *polygon,
                              gpc_tristrip    *tristrip);

void gpc_free_polygon        (gpc_polygon     *polygon);

void gpc_free_tristrip       (gpc_tristrip    *tristrip);

#endif // !defined(AFX_POLYGON_H__6B661CC4_2E36_4F4F_8F4D_872D408901D6__INCLUDED_)
