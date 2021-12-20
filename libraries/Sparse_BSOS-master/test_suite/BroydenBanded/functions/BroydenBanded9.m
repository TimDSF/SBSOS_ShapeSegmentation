%% This file generates the Broyden Banded function in 9 variables 
 sf =10^-5; 
 scale = false; 
 
%% generate polynomial F  nvar = 9.0;
  F = sparse(233,9);
  F(1,1) = 6.0;
  F(2,1) = 4.0;
  F(3,1) = 3.0;
  F(3,2) = 2.0;
  F(4,1) = 3.0;
  F(4,2) = 1.0;
  F(5,1) = 3.0;
  F(6,1) = 2.0;
  F(6,2) = 3.0;
  F(7,1) = 2.0;
  F(7,2) = 2.0;
  F(8,1) = 2.0;
  F(8,2) = 1.0;
  F(9,1) = 2.0;
  F(9,3) = 3.0;
  F(10,1) = 2.0;
  F(10,3) = 2.0;
  F(11,1) = 2.0;
  F(11,3) = 1.0;
  F(12,1) = 2.0;
  F(12,4) = 3.0;
  F(13,1) = 2.0;
  F(13,4) = 2.0;
  F(14,1) = 2.0;
  F(14,4) = 1.0;
  F(15,1) = 2.0;
  F(15,5) = 3.0;
  F(16,1) = 2.0;
  F(16,5) = 2.0;
  F(17,1) = 2.0;
  F(17,6) = 3.0;
  F(18,1) = 2.0;
  F(18,6) = 2.0;
  F(19,1) = 2.0;
  F(19,6) = 1.0;
  F(20,1) = 2.0;
  F(20,7) = 2.0;
  F(21,1) = 2.0;
  F(21,7) = 1.0;
  F(22,1) = 2.0;
  F(23,1) = 1.0;
  F(23,2) = 3.0;
  F(24,1) = 1.0;
  F(24,2) = 2.0;
  F(25,1) = 1.0;
  F(25,3) = 3.0;
  F(26,1) = 1.0;
  F(26,3) = 2.0;
  F(27,1) = 1.0;
  F(27,3) = 1.0;
  F(28,1) = 1.0;
  F(28,4) = 3.0;
  F(29,1) = 1.0;
  F(29,4) = 2.0;
  F(30,1) = 1.0;
  F(30,4) = 1.0;
  F(31,1) = 1.0;
  F(31,5) = 3.0;
  F(32,1) = 1.0;
  F(32,5) = 2.0;
  F(33,1) = 1.0;
  F(33,6) = 3.0;
  F(34,1) = 1.0;
  F(34,6) = 2.0;
  F(35,1) = 1.0;
  F(35,6) = 1.0;
  F(36,1) = 1.0;
  F(36,7) = 2.0;
  F(37,1) = 1.0;
  F(37,7) = 1.0;
  F(38,1) = 1.0;
  F(39,2) = 6.0;
  F(40,2) = 4.0;
  F(41,2) = 3.0;
  F(41,3) = 2.0;
  F(42,2) = 3.0;
  F(42,3) = 1.0;
  F(43,2) = 3.0;
  F(44,2) = 2.0;
  F(44,3) = 3.0;
  F(45,2) = 2.0;
  F(45,3) = 2.0;
  F(46,2) = 2.0;
  F(46,3) = 1.0;
  F(47,2) = 2.0;
  F(47,4) = 3.0;
  F(48,2) = 2.0;
  F(48,4) = 2.0;
  F(49,2) = 2.0;
  F(49,4) = 1.0;
  F(50,2) = 2.0;
  F(50,5) = 3.0;
  F(51,2) = 2.0;
  F(51,5) = 2.0;
  F(52,2) = 2.0;
  F(52,5) = 1.0;
  F(53,2) = 2.0;
  F(53,6) = 3.0;
  F(54,2) = 2.0;
  F(54,6) = 2.0;
  F(55,2) = 2.0;
  F(55,7) = 3.0;
  F(56,2) = 2.0;
  F(56,7) = 2.0;
  F(57,2) = 2.0;
  F(57,7) = 1.0;
  F(58,2) = 2.0;
  F(58,8) = 2.0;
  F(59,2) = 2.0;
  F(59,8) = 1.0;
  F(60,2) = 2.0;
  F(61,2) = 1.0;
  F(61,3) = 3.0;
  F(62,2) = 1.0;
  F(62,3) = 2.0;
  F(63,2) = 1.0;
  F(63,4) = 3.0;
  F(64,2) = 1.0;
  F(64,4) = 2.0;
  F(65,2) = 1.0;
  F(65,4) = 1.0;
  F(66,2) = 1.0;
  F(66,5) = 3.0;
  F(67,2) = 1.0;
  F(67,5) = 2.0;
  F(68,2) = 1.0;
  F(68,5) = 1.0;
  F(69,2) = 1.0;
  F(69,6) = 3.0;
  F(70,2) = 1.0;
  F(70,6) = 2.0;
  F(71,2) = 1.0;
  F(71,7) = 3.0;
  F(72,2) = 1.0;
  F(72,7) = 2.0;
  F(73,2) = 1.0;
  F(73,7) = 1.0;
  F(74,2) = 1.0;
  F(74,8) = 2.0;
  F(75,2) = 1.0;
  F(75,8) = 1.0;
  F(76,2) = 1.0;
  F(77,3) = 6.0;
  F(78,3) = 4.0;
  F(79,3) = 3.0;
  F(79,4) = 2.0;
  F(80,3) = 3.0;
  F(80,4) = 1.0;
  F(81,3) = 3.0;
  F(82,3) = 2.0;
  F(82,4) = 3.0;
  F(83,3) = 2.0;
  F(83,4) = 2.0;
  F(84,3) = 2.0;
  F(84,4) = 1.0;
  F(85,3) = 2.0;
  F(85,5) = 3.0;
  F(86,3) = 2.0;
  F(86,5) = 2.0;
  F(87,3) = 2.0;
  F(87,5) = 1.0;
  F(88,3) = 2.0;
  F(88,6) = 3.0;
  F(89,3) = 2.0;
  F(89,6) = 2.0;
  F(90,3) = 2.0;
  F(90,6) = 1.0;
  F(91,3) = 2.0;
  F(91,7) = 3.0;
  F(92,3) = 2.0;
  F(92,7) = 2.0;
  F(93,3) = 2.0;
  F(93,8) = 3.0;
  F(94,3) = 2.0;
  F(94,8) = 2.0;
  F(95,3) = 2.0;
  F(95,8) = 1.0;
  F(96,3) = 2.0;
  F(96,9) = 2.0;
  F(97,3) = 2.0;
  F(97,9) = 1.0;
  F(98,3) = 2.0;
  F(99,3) = 1.0;
  F(99,4) = 3.0;
  F(100,3) = 1.0;
  F(100,4) = 2.0;
  F(101,3) = 1.0;
  F(101,5) = 3.0;
  F(102,3) = 1.0;
  F(102,5) = 2.0;
  F(103,3) = 1.0;
  F(103,5) = 1.0;
  F(104,3) = 1.0;
  F(104,6) = 3.0;
  F(105,3) = 1.0;
  F(105,6) = 2.0;
  F(106,3) = 1.0;
  F(106,6) = 1.0;
  F(107,3) = 1.0;
  F(107,7) = 3.0;
  F(108,3) = 1.0;
  F(108,7) = 2.0;
  F(109,3) = 1.0;
  F(109,8) = 3.0;
  F(110,3) = 1.0;
  F(110,8) = 2.0;
  F(111,3) = 1.0;
  F(111,8) = 1.0;
  F(112,3) = 1.0;
  F(112,9) = 2.0;
  F(113,3) = 1.0;
  F(113,9) = 1.0;
  F(114,3) = 1.0;
  F(115,4) = 6.0;
  F(116,4) = 4.0;
  F(117,4) = 3.0;
  F(117,5) = 2.0;
  F(118,4) = 3.0;
  F(118,5) = 1.0;
  F(119,4) = 3.0;
  F(120,4) = 2.0;
  F(120,5) = 3.0;
  F(121,4) = 2.0;
  F(121,5) = 2.0;
  F(122,4) = 2.0;
  F(122,5) = 1.0;
  F(123,4) = 2.0;
  F(123,6) = 3.0;
  F(124,4) = 2.0;
  F(124,6) = 2.0;
  F(125,4) = 2.0;
  F(125,6) = 1.0;
  F(126,4) = 2.0;
  F(126,7) = 3.0;
  F(127,4) = 2.0;
  F(127,7) = 2.0;
  F(128,4) = 2.0;
  F(128,7) = 1.0;
  F(129,4) = 2.0;
  F(129,8) = 3.0;
  F(130,4) = 2.0;
  F(130,8) = 2.0;
  F(131,4) = 2.0;
  F(131,9) = 3.0;
  F(132,4) = 2.0;
  F(132,9) = 2.0;
  F(133,4) = 2.0;
  F(133,9) = 1.0;
  F(134,4) = 2.0;
  F(135,4) = 1.0;
  F(135,5) = 3.0;
  F(136,4) = 1.0;
  F(136,5) = 2.0;
  F(137,4) = 1.0;
  F(137,6) = 3.0;
  F(138,4) = 1.0;
  F(138,6) = 2.0;
  F(139,4) = 1.0;
  F(139,6) = 1.0;
  F(140,4) = 1.0;
  F(140,7) = 3.0;
  F(141,4) = 1.0;
  F(141,7) = 2.0;
  F(142,4) = 1.0;
  F(142,7) = 1.0;
  F(143,4) = 1.0;
  F(143,8) = 3.0;
  F(144,4) = 1.0;
  F(144,8) = 2.0;
  F(145,4) = 1.0;
  F(145,9) = 3.0;
  F(146,4) = 1.0;
  F(146,9) = 2.0;
  F(147,4) = 1.0;
  F(147,9) = 1.0;
  F(148,4) = 1.0;
  F(149,5) = 6.0;
  F(150,5) = 4.0;
  F(151,5) = 3.0;
  F(151,6) = 2.0;
  F(152,5) = 3.0;
  F(152,6) = 1.0;
  F(153,5) = 3.0;
  F(154,5) = 2.0;
  F(154,6) = 3.0;
  F(155,5) = 2.0;
  F(155,6) = 2.0;
  F(156,5) = 2.0;
  F(156,6) = 1.0;
  F(157,5) = 2.0;
  F(157,7) = 3.0;
  F(158,5) = 2.0;
  F(158,7) = 2.0;
  F(159,5) = 2.0;
  F(159,7) = 1.0;
  F(160,5) = 2.0;
  F(160,8) = 3.0;
  F(161,5) = 2.0;
  F(161,8) = 2.0;
  F(162,5) = 2.0;
  F(162,9) = 3.0;
  F(163,5) = 2.0;
  F(163,9) = 2.0;
  F(164,5) = 2.0;
  F(164,9) = 1.0;
  F(165,5) = 2.0;
  F(166,5) = 1.0;
  F(166,6) = 3.0;
  F(167,5) = 1.0;
  F(167,6) = 2.0;
  F(168,5) = 1.0;
  F(168,6) = 1.0;
  F(169,5) = 1.0;
  F(169,7) = 3.0;
  F(170,5) = 1.0;
  F(170,7) = 2.0;
  F(171,5) = 1.0;
  F(171,7) = 1.0;
  F(172,5) = 1.0;
  F(172,8) = 3.0;
  F(173,5) = 1.0;
  F(173,8) = 2.0;
  F(174,5) = 1.0;
  F(174,9) = 3.0;
  F(175,5) = 1.0;
  F(175,9) = 2.0;
  F(176,5) = 1.0;
  F(176,9) = 1.0;
  F(177,5) = 1.0;
  F(178,6) = 6.0;
  F(179,6) = 4.0;
  F(180,6) = 3.0;
  F(180,7) = 2.0;
  F(181,6) = 3.0;
  F(181,7) = 1.0;
  F(182,6) = 3.0;
  F(183,6) = 2.0;
  F(183,7) = 3.0;
  F(184,6) = 2.0;
  F(184,7) = 2.0;
  F(185,6) = 2.0;
  F(185,8) = 3.0;
  F(186,6) = 2.0;
  F(186,8) = 2.0;
  F(187,6) = 2.0;
  F(187,9) = 3.0;
  F(188,6) = 2.0;
  F(188,9) = 2.0;
  F(189,6) = 2.0;
  F(189,9) = 1.0;
  F(190,6) = 1.0;
  F(190,7) = 3.0;
  F(191,6) = 1.0;
  F(191,7) = 1.0;
  F(192,6) = 1.0;
  F(192,8) = 3.0;
  F(193,6) = 1.0;
  F(193,8) = 2.0;
  F(194,6) = 1.0;
  F(194,9) = 3.0;
  F(195,6) = 1.0;
  F(195,9) = 2.0;
  F(196,6) = 1.0;
  F(196,9) = 1.0;
  F(197,6) = 1.0;
  F(198,7) = 6.0;
  F(199,7) = 4.0;
  F(200,7) = 3.0;
  F(200,8) = 2.0;
  F(201,7) = 3.0;
  F(201,8) = 1.0;
  F(202,7) = 3.0;
  F(203,7) = 2.0;
  F(203,8) = 3.0;
  F(204,7) = 2.0;
  F(204,8) = 2.0;
  F(205,7) = 2.0;
  F(205,8) = 1.0;
  F(206,7) = 2.0;
  F(206,9) = 3.0;
  F(207,7) = 2.0;
  F(207,9) = 2.0;
  F(208,7) = 2.0;
  F(208,9) = 1.0;
  F(209,7) = 2.0;
  F(210,7) = 1.0;
  F(210,8) = 3.0;
  F(211,7) = 1.0;
  F(211,8) = 2.0;
  F(212,7) = 1.0;
  F(212,8) = 1.0;
  F(213,7) = 1.0;
  F(213,9) = 3.0;
  F(214,7) = 1.0;
  F(214,9) = 2.0;
  F(215,7) = 1.0;
  F(215,9) = 1.0;
  F(216,7) = 1.0;
  F(217,8) = 6.0;
  F(218,8) = 4.0;
  F(219,8) = 3.0;
  F(219,9) = 2.0;
  F(220,8) = 3.0;
  F(220,9) = 1.0;
  F(221,8) = 3.0;
  F(222,8) = 2.0;
  F(222,9) = 3.0;
  F(223,8) = 2.0;
  F(223,9) = 1.0;
  F(224,8) = 2.0;
  F(225,8) = 1.0;
  F(225,9) = 3.0;
  F(226,8) = 1.0;
  F(226,9) = 2.0;
  F(227,8) = 1.0;
  F(227,9) = 1.0;
  F(228,9) = 6.0;
  F(229,9) = 4.0;
  F(230,9) = 3.0;
  F(231,9) = 2.0;
  F(232,9) = 1.0;
  %Coefficients:
  F(1,10) = 1.0e2;
  F(2,10) = 4.5e1;
  F(3,10) = -2.0e1;
  F(4,10) = -2.0e1;
  F(5,10) = 3.0e1;
  F(6,10) = -2.0e1;
  F(7,10) = 8.0;
  F(8,10) = 4.0;
  F(9,10) = -2.0e1;
  F(10,10) = 8.0;
  F(11,10) = 4.0;
  F(12,10) = -2.0e1;
  F(13,10) = 6.0;
  F(14,10) = 2.0;
  F(15,10) = -2.0e1;
  F(16,10) = 4.0;
  F(17,10) = -2.0e1;
  F(18,10) = 2.0;
  F(19,10) = -2.0;
  F(20,10) = 2.0;
  F(21,10) = 2.0;
  F(22,10) = -1.0;
  F(23,10) = -2.0e1;
  F(24,10) = 4.0;
  F(25,10) = -2.0e1;
  F(26,10) = 8.0;
  F(27,10) = 4.0;
  F(28,10) = -2.0e1;
  F(29,10) = 6.0;
  F(30,10) = 2.0;
  F(31,10) = -2.0e1;
  F(32,10) = 4.0;
  F(33,10) = -2.0e1;
  F(34,10) = 2.0;
  F(35,10) = -2.0;
  F(36,10) = 2.0;
  F(37,10) = 2.0;
  F(38,10) = -6.0;
  F(39,10) = 1.0e2;
  F(40,10) = 4.6e1;
  F(41,10) = -2.0e1;
  F(42,10) = -2.0e1;
  F(43,10) = 3.2e1;
  F(44,10) = -2.0e1;
  F(45,10) = 8.0;
  F(46,10) = 4.0;
  F(47,10) = -2.0e1;
  F(48,10) = 8.0;
  F(49,10) = 4.0;
  F(50,10) = -2.0e1;
  F(51,10) = 6.0;
  F(52,10) = 2.0;
  F(53,10) = -2.0e1;
  F(54,10) = 4.0;
  F(55,10) = -2.0e1;
  F(56,10) = 2.0;
  F(57,10) = -2.0;
  F(58,10) = 2.0;
  F(59,10) = 2.0;
  F(60,10) = -2.0;
  F(61,10) = -2.0e1;
  F(62,10) = 4.0;
  F(63,10) = -2.0e1;
  F(64,10) = 8.0;
  F(65,10) = 4.0;
  F(66,10) = -2.0e1;
  F(67,10) = 6.0;
  F(68,10) = 2.0;
  F(69,10) = -2.0e1;
  F(70,10) = 4.0;
  F(71,10) = -2.0e1;
  F(72,10) = 2.0;
  F(73,10) = -2.0;
  F(74,10) = 2.0;
  F(75,10) = 2.0;
  F(76,10) = -8.0;
  F(77,10) = 1.0e2;
  F(78,10) = 4.6e1;
  F(79,10) = -2.0e1;
  F(80,10) = -2.0e1;
  F(81,10) = 3.2e1;
  F(82,10) = -2.0e1;
  F(83,10) = 8.0;
  F(84,10) = 4.0;
  F(85,10) = -2.0e1;
  F(86,10) = 8.0;
  F(87,10) = 4.0;
  F(88,10) = -2.0e1;
  F(89,10) = 6.0;
  F(90,10) = 2.0;
  F(91,10) = -2.0e1;
  F(92,10) = 4.0;
  F(93,10) = -2.0e1;
  F(94,10) = 2.0;
  F(95,10) = -2.0;
  F(96,10) = 2.0;
  F(97,10) = 2.0;
  F(98,10) = -2.0;
  F(99,10) = -2.0e1;
  F(100,10) = 4.0;
  F(101,10) = -2.0e1;
  F(102,10) = 8.0;
  F(103,10) = 4.0;
  F(104,10) = -2.0e1;
  F(105,10) = 6.0;
  F(106,10) = 2.0;
  F(107,10) = -2.0e1;
  F(108,10) = 4.0;
  F(109,10) = -2.0e1;
  F(110,10) = 2.0;
  F(111,10) = -2.0;
  F(112,10) = 2.0;
  F(113,10) = 2.0;
  F(114,10) = -8.0;
  F(115,10) = 1.0e2;
  F(116,10) = 4.6e1;
  F(117,10) = -2.0e1;
  F(118,10) = -2.0e1;
  F(119,10) = 3.2e1;
  F(120,10) = -2.0e1;
  F(121,10) = 8.0;
  F(122,10) = 4.0;
  F(123,10) = -2.0e1;
  F(124,10) = 8.0;
  F(125,10) = 4.0;
  F(126,10) = -2.0e1;
  F(127,10) = 6.0;
  F(128,10) = 2.0;
  F(129,10) = -2.0e1;
  F(130,10) = 4.0;
  F(131,10) = -2.0e1;
  F(132,10) = 2.0;
  F(133,10) = -2.0;
  F(134,10) = -2.0;
  F(135,10) = -2.0e1;
  F(136,10) = 4.0;
  F(137,10) = -2.0e1;
  F(138,10) = 8.0;
  F(139,10) = 4.0;
  F(140,10) = -2.0e1;
  F(141,10) = 6.0;
  F(142,10) = 2.0;
  F(143,10) = -2.0e1;
  F(144,10) = 4.0;
  F(145,10) = -2.0e1;
  F(146,10) = 2.0;
  F(147,10) = -2.0;
  F(148,10) = -8.0;
  F(149,10) = 1.0e2;
  F(150,10) = 4.5e1;
  F(151,10) = -2.0e1;
  F(152,10) = -2.0e1;
  F(153,10) = 3.0e1;
  F(154,10) = -2.0e1;
  F(155,10) = 6.0;
  F(156,10) = 2.0;
  F(157,10) = -2.0e1;
  F(158,10) = 6.0;
  F(159,10) = 2.0;
  F(160,10) = -2.0e1;
  F(161,10) = 4.0;
  F(162,10) = -2.0e1;
  F(163,10) = 2.0;
  F(164,10) = -2.0;
  F(165,10) = -1.0;
  F(166,10) = -2.0e1;
  F(167,10) = 2.0;
  F(168,10) = -2.0;
  F(169,10) = -2.0e1;
  F(170,10) = 6.0;
  F(171,10) = 2.0;
  F(172,10) = -2.0e1;
  F(173,10) = 4.0;
  F(174,10) = -2.0e1;
  F(175,10) = 2.0;
  F(176,10) = -2.0;
  F(177,10) = -6.0;
  F(178,10) = 1.0e2;
  F(179,10) = 4.4e1;
  F(180,10) = -2.0e1;
  F(181,10) = -2.0e1;
  F(182,10) = 2.8e1;
  F(183,10) = -2.0e1;
  F(184,10) = 4.0;
  F(185,10) = -2.0e1;
  F(186,10) = 4.0;
  F(187,10) = -2.0e1;
  F(188,10) = 2.0;
  F(189,10) = -2.0;
  F(190,10) = -2.0e1;
  F(191,10) = -4.0;
  F(192,10) = -2.0e1;
  F(193,10) = 4.0;
  F(194,10) = -2.0e1;
  F(195,10) = 2.0;
  F(196,10) = -2.0;
  F(197,10) = -4.0;
  F(198,10) = 1.0e2;
  F(199,10) = 4.3e1;
  F(200,10) = -2.0e1;
  F(201,10) = -2.0e1;
  F(202,10) = 2.6e1;
  F(203,10) = -2.0e1;
  F(204,10) = 2.0;
  F(205,10) = -2.0;
  F(206,10) = -2.0e1;
  F(207,10) = 2.0;
  F(208,10) = -2.0;
  F(209,10) = 1.0;
  F(210,10) = -2.0e1;
  F(211,10) = -2.0;
  F(212,10) = -6.0;
  F(213,10) = -2.0e1;
  F(214,10) = 2.0;
  F(215,10) = -2.0;
  F(216,10) = -2.0;
  F(217,10) = 1.0e2;
  F(218,10) = 4.2e1;
  F(219,10) = -2.0e1;
  F(220,10) = -2.0e1;
  F(221,10) = 2.4e1;
  F(222,10) = -2.0e1;
  F(223,10) = -4.0;
  F(224,10) = 2.0;
  F(225,10) = -2.0e1;
  F(226,10) = -4.0;
  F(227,10) = -8.0;
  F(228,10) = 1.0e2;
  F(229,10) = 4.1e1;
  F(230,10) = 2.2e1;
  F(231,10) = 3.0;
  F(232,10) = 2.0;
  F(233,10) = 9.0;

 if scale 
 F(:,9+1) = F(:,9+1)*sf; 
 end
 G ={};
 
I={[1, 2, 3, 4, 5, 6, 7], [2, 3, 4, 5, 6, 7, 8], [3, 4, 5, 6, 7, 8, 9]};
 
J={[], [], []};

pop.F = F; pop.I = I; pop.G = G; pop.J = J; pop.n = size(F,2)-1;