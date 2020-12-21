/******************************************************************************
 * Project:  libngstore
 * Purpose:  NextGIS store and visualisation support library
 * Author: Dmitry Baryshnikov, dmitry.baryshnikov@nextgis.com
 ******************************************************************************
 *   Copyright (c) 2017-2020 NextGIS, <info@nextgis.com>
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU Lesser General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#ifndef NGUVERSION_H
#define NGUVERSION_H

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#define NGU_VERSION_MAJOR 0
#define NGU_VERSION_MINOR 2
#define NGU_VERSION_REV   1
#define NGU_VERSION  STR(NGU_VERSION_MAJOR) "." STR(NGU_VERSION_MINOR) "." \
    STR(NGU_VERSION_REV)

#define NGU_COMPUTE_VERSION(maj,min,rev) ((maj)*10000+(min)*100+rev) // maj - any, min < 99, rev < 99
#define NGU_VERSION_NUM NGU_COMPUTE_VERSION(NGU_VERSION_MAJOR,NGU_VERSION_MINOR,NGU_VERSION_REV)

/*  check if the current version is at least major.minor.revision */
#define CHECK_VERSION(major,minor,rev) \
    (NGU_VERSION_MAJOR > (major) || \
    (NGU_VERSION_MAJOR == (major) && NGU_VERSION_MINOR > (minor)) || \
    (NGU_VERSION_MAJOR == (major) && NGU_VERSION_MINOR == (minor) && \
     NGU_VERSION_REV >= (release)))

#endif // NGUVERSION_H
