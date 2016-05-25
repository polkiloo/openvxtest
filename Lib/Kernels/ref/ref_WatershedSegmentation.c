/*
File: ref_watershedSegmentation.c
Iplementation of watershed segmentation algorithm
Author: Levitskiy Mikhail
Date: 27 April 2016
*/
#include "../ref.h"
#include <assert.h>
#define ABS(a)  ( (a)>=0 ? (a) : -(a) )
// Labels for pixels
const vx_int8 IN_QUEUE = -2; // Pixel visited
const vx_int8 WSHED = -1;    // Pixel belongs to watershed
const vx_uint16 NQ = 256; // possible bit values = 2^8
const vx_uint8 NCH = 3; // number of channels of the input image
//Vx Watershed Nodes
typedef struct VxWSNode
{
   vx_uint32 next;
   vx_int32 mask_ofs;
   vx_int32 img_ofs;
}
VxWSNode;
// Queue for WSNodes
typedef struct VxWSQueue
{
   vx_uint32 first;
   vx_uint32 last;
}
VxWSQueue;
static vx_uint32   allocVxWSNodes(VxWSNode** storage, vx_uint32* size);
static bool        checkImagesParametrs(const vx_image* src_image, const vx_image* markers);
//Watershed segmentation
vx_status          ref_WatershedSegmentation(const vx_image src_image, vx_image markers)
{
   //check images param  and continue  if it`s valid
   if (checkImagesParametrs(&src_image, &markers))
   {
      return VX_ERROR_INVALID_PARAMETERS;
   }
   //size image
   vx_uint32 img_width = src_image->width;
   vx_uint32 img_height = src_image->height;
   // storage size
   vx_uint32 storageSize = 0;
   //array of every created node
   struct VxWSNode* storage = (VxWSNode*)malloc(sizeof (VxWSNode));
   vx_uint32 free_node = 0, node;
   // Priority queue of queues of nodes
   // from high priority (0) to low priority (255)
   VxWSQueue* q = (VxWSQueue*)calloc(NQ, sizeof(VxWSQueue));
   // Non-empty queue with highest priority
   vx_uint32 active_queue;
   vx_uint32 i, j;
   // Color differences
   vx_uint8 dr, dg, db;
   vx_uint32 subs_tab[513];
   for (i = 0; i < 256; i++)
      subs_tab[i] = 0;
   for (i = 256; i <= 512; i++)
      subs_tab[i] = i - 256;
   // Current pixel in input image
   const vx_uint8* img = (vx_uint8*)(src_image->data);
   // Step size to next row in input image
   vx_int32 istep = (vx_int32)(img_width * NCH);
   // Current pixel in mask image
   vx_int32* mask = (vx_int32*)(markers->data);
   // Step size to next row in mask image
   vx_int32 mstep = (vx_int32)(img_width);
   /////////////////#INLINE FUNCTION///////////////////////////////////////
#define ws_max(a,b) ((b) + subs_tab[(a)-(b)+NQ])
#define ws_min(a,b) ((a) - subs_tab[(a)-(b)+NQ])
   // Create a new node with offsets mofs and iofs in queue idx 
#define ws_push(idx,mofs,iofs)                          \
   {                                                    \
   if (!free_node)                                      \
   free_node = allocVxWSNodes(&storage, &storageSize);  \
   node = free_node;                                    \
   free_node = storage[free_node].next;                 \
   storage[node].next = 0;                              \
   storage[node].mask_ofs = mofs;                       \
   storage[node].img_ofs = iofs;                        \
   if (q[idx].last)                                     \
   storage[q[idx].last].next = node;                    \
   else                                                 \
   q[idx].first = node;                                 \
   q[idx].last = node;                                  \
   }
   // Get next node from queue idx
#define ws_pop(idx,mofs,iofs)                           \
   {                                                    \
   node = q[idx].first;                                 \
   q[idx].first = storage[node].next;                   \
   if (!storage[node].next)                             \
   q[idx].last = 0;                                     \
   storage[node].next = free_node;                      \
   free_node = node;                                    \
   mofs = storage[node].mask_ofs;                       \
   iofs = storage[node].img_ofs;                        \
   }  \
   // Get highest absolute channel difference in diff
#define c_diff(ptr1,ptr2,diff)                          \
   {                                                    \
   dr = ABS((ptr1)[0] - (ptr2)[0]);                     \
   dg = ABS((ptr1)[1] - (ptr2)[1]);                     \
   db = ABS((ptr1)[2] - (ptr2)[2]);                     \
   diff = ws_max(db, dg);                               \
   diff = ws_max(diff, dr);                             \
   }
   ////////////////////////////////////////////////////////////////////////
   // draw a pixel-wide border of dummy "watershed" (i.e. boundary) pixels
   for (j = 0; j < img_width; j++)
      mask[j] = mask[j + mstep*(img_height - 1)] = WSHED;
   // initial phase: put all the neighbor pixels of each marker to the ordered queue -
   // determine the initial boundaries of the basins
   for (i = 1; i < img_height - 1; i++)
   {
      img += istep; mask += mstep;
      mask[0] = mask[img_width - 1] = WSHED; //boundary pixels
      for (j = 1; j <img_width - 1; j++)
      {
         vx_int32* m = mask + j;
         if (m[0] < 0) m[0] = 0;
         if (m[0] == 0 && (m[-1] > 0 || m[1] > 0 || m[-mstep] > 0 || m[mstep] > 0))
         {
            // Find smallest difference to adjacent markers
            const vx_uint8* ptr = img + j * NCH;
            vx_uint32 idx = 256, t;
            if (m[-1] > 0)
               c_diff(ptr, ptr - NCH, idx);
            if (m[1] > 0)
            {
               c_diff(ptr, ptr + NCH, t);
               idx = ws_min(idx, t);
            }
            if (m[-mstep] > 0)
            {
               c_diff(ptr, ptr - istep, t);
               idx = ws_min(idx, t);
            }
            if (m[mstep] > 0)
            {
               c_diff(ptr, ptr + istep, t);
               idx = ws_min(idx, t);
            }
            // Add to according queue
            assert(0 <= idx && idx <= 255);
            ws_push(idx, i*mstep + j, i*istep + j * NCH);
            m[0] = IN_QUEUE;
         }
      }
   }
   for (i = 0; i < NQ; i++)
   if (q[i].first)
      break;
   // if there is no markers, exit immediately
   if (i == NQ)
      return VX_ERROR_INVALID_PARAMETERS;
   active_queue = i;
   img = (vx_uint8*)(src_image->data);
   mask = (vx_int32*)(markers->data);
   // recursively fill the basins
   for (;;)
   {
      vx_int32 mofs, iofs;
      vx_int32 lab = 0, t;
      vx_int32* m;
      const vx_uint8* ptr;
      // Get non-empty queue with highest priority
      // Exit condition: empty priority queue
      if (q[active_queue].first == 0)
      {
         for (i = active_queue + 1; i < NQ; i++)
         if (q[i].first)
            break;
         if (i == NQ)
            break;
         active_queue = i;
      }
      // Get next node
      ws_pop(active_queue, mofs, iofs);
      // Calculate pointer to current pixel in input and marker image
      m = mask + mofs;
      ptr = img + iofs;
      // Check surrounding pixels for labels
      // to determine label for current pixel
      t = m[-1]; // Left
      if (t > 0) lab = t;
      t = m[1]; // Right
      if (t > 0)
      {
         if (lab == 0) lab = t;
         else if (t != lab) lab = WSHED;
      }
      t = m[-mstep]; // Top
      if (t > 0)
      {
         if (lab == 0) lab = t;
         else if (t != lab) lab = WSHED;
      }
      t = m[mstep]; // Bottom
      if (t > 0)
      {
         if (lab == 0) lab = t;
         else if (t != lab) lab = WSHED;
      }
      // Set label to current pixel in marker image
      //assert(wslab != 0);
      m[0] = lab;
      if (lab == WSHED)
         continue;
      // Add adjacent, unlabeled pixels to corresponding queue
      if (m[-1] == 0)
      {
         c_diff(ptr, ptr - NCH, t);
         ws_push(t, mofs - 1, iofs - NCH);
         active_queue = ws_min(active_queue, t);
         m[-1] = IN_QUEUE;
      }
      if (m[1] == 0)
      {
         c_diff(ptr, ptr + NCH, t);
         ws_push(t, mofs + 1, iofs + NCH);
         active_queue = ws_min(active_queue, t);
         m[1] = IN_QUEUE;
      }
      if (m[-mstep] == 0)
      {
         c_diff(ptr, ptr - istep, t);
         ws_push(t, mofs - mstep, iofs - istep);
         active_queue = ws_min(active_queue, t);
         m[-mstep] = IN_QUEUE;
      }
      if (m[mstep] == 0)
      {
         c_diff(ptr, ptr + istep, t);
         ws_push(t, mofs + mstep, iofs + istep);
         active_queue = ws_min(active_queue, t);
         m[mstep] = IN_QUEUE;
      }
   }
   free(storage);
   free(q);
   return VX_SUCCESS;
}
//Allocate nodes in storage of VxWSNode
static vx_uint32   allocVxWSNodes(VxWSNode** storage, vx_uint32* size)
{
   vx_uint32 sz = *size;
   vx_uint32 newsz = (128 < (sz * 3 / 2) ? (sz * 3 / 2) : 128);
   *storage = (VxWSNode*)realloc(*storage, (newsz)*sizeof(VxWSNode));
   *size = newsz;
   if (sz == 0)
   {
      (*storage)[0].next = 0;
      sz = 1;
   }
   for (vx_uint32 i = sz; i < newsz - 1; i++)
   {
      (*storage)[i].next = i + 1;
   }
   (*storage)[newsz - 1].next = 0;
   return sz;
}
//Check input images format
static bool        checkImagesParametrs(const vx_image* src_image, const vx_image* markers)
{
   return (
                ((*src_image)->image_type != VX_DF_IMAGE_RGB) 
             || ((*markers)->image_type   != VX_DF_IMAGE_S32) 
             || ((*src_image)->height     != (*markers)->height) 
             || ((*src_image)->width      != (*markers)->width)
          );
}
